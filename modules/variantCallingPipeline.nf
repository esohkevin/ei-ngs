def getBamFileSet() {
    return channel.fromFilePairs( params.output_dir + '/bam/' + "*.{bam,bam.bai}", size: 2, flat: true )
                  .ifEmpty { error "\nERROR: Could not locate a file! \n" }
                  .map { bamName, bamFile, bamIndex -> tuple(bamName, bamFile, bamIndex) }
}

def getCramFileSet() {
    return channel.fromFilePairs( params.output_dir + '/cram/' + "*.{cram,cram.crai}", size: 2, flat: true )
                  .ifEmpty { error "\nERROR: Could not locate a file! \n" }
                  .map { bamName, bamFile, bamIndex -> tuple(bamName, bamFile, bamIndex) }
}

def getAlignmentFileSet() {
    return channel.fromFilePairs( [ params.alignment_dir + "/*.{bam,bam.bai}", params.alignment_dir + "*.{cram,cram.crai}" ] , size: 2, flat: true )
                  .ifEmpty { error "\nERROR: Could not locate a file! \n" }
                  .map { bamName, bamFile, bamIndex -> tuple(bamName, bamFile, bamIndex) }
}

def getGvcfFiles() {
    return channel.fromPath( params.gvcf_dir + "/*.{gvcf,g.vcf}.{gz,gz.tbi}" )
                  .flatten()
}

def getGenomicsdbWorkspaces() {
    return channel.fromPath( params.genomicsdb_workspace_dir + "/*", type: 'dir' )
                  .flatten()
}


def getGenomicInterval(gvcfList) {
    if(params.interval == "NULL") {
        genomicInterval = getVcfGenomicIntervals(gvcfList).flatten()
    }
    else {
        genomicInterval = getGenomicIntervalList().flatten()
    }
}


process getGvcfList() {
    tag "creating GVCF list..."
    input:
        path(gvcfList)
    output:
        path("gvcf.list")
    script:
        """
        readlink *.gz | awk '{print "-V",\$1}' > gvcf.list

        #for file in ${gvcfList}; do
        #    if [ \${file##*.} != "tbi" ]; then
        #        echo "-V \$(readlink \${file})";
        #    fi
        #done > gvcf.list        
        """
}

// TO RUN HAPLOTYPE CALLER PER INTERVAL
process getAlignmentGenomicIntervals() {
    tag "extracting invervals from alignment file (${bamName})..."
    label 'samtools'
    label 'smallMemory'
    cache 'lenient'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        path("*.bed")
    script:
        """
        ######################################################
        # make non-overlapping intervals of 5,000,000 bp for #
        # chromosomes larger than 5,000,000. Otherwise, just #
        # print out chromosome length as only interval.      #
        ######################################################

        ### BED INTERVAL
        samtools \
            view \
            -H ${bamFile} | \
        grep 'SQ' | \
        cut -f2-3 | \
        sed 's/[SL]N://g' | \
        awk '{ if(\$2<=5000000){print \$1,"0",\$2} else{ for(i=0; i<=\$2; i+=5000000) { if(i+4999999<\$2) {print \$1,i,i+4999999} else{print \$1,i,\$2} } } }' \
        > .interval_list

        while read interval; do
            echo \$interval > \$(echo \${interval} | sed 's/[:*]/_/g' | sed 's/ /_/g').bed
        done < .interval_list

        ### GATK LIST INTERVAL
        #samtools \
        #    view \
        #    -H ${bamFile} | \
        #grep 'SQ' | \
        #cut -f2-3 | \
        #sed 's/[SL]N://g' | \
        #        awk '{ if(\$2<=5000000){print \$1":0-"\$2} else{ for(i=0; i<=\$2; i+=5000000) { if(i+4999999<\$2) {print \$1":"i"-"i+4999999} else{print \$1":"i"-"\$2}} } }' \
        #for interval in \$(cat .interval_list); do
        #    echo \$interval > \$(echo \${interval} | sed 's/[:*]/_/g').bed
        #done
        """
}

// TO COMBINE PER SAMPLE VARIANT CALLS AND 
// RUN JOINT VARIANT CALLING PER INTERVAL
process getVcfGenomicIntervals() {
    tag "extracting invervals from GVCF file..."
    label 'variantCaller'
    cache 'lenient'
    input:
        path gvcfList
    output:
        path "*"
    script:
        """
        gvcfile=\$(readlink ${gvcfList})
        gvcf=\$(awk '{print \$2}' \${gvcfile} | head -1)

        zgrep '##contig' \${gvcf} | \
            sed 's/[=,>]/\t/g' | \
            cut -f3,5 | \
            grep -v '^HLA' | \
        awk '{ if(\$2<=5000000){print \$1,"0",\$2} else{ for(i=0; i<=\$2; i+=5000000) { if(i+4999999<\$2) {print \$1,i,i+4999999} else{print \$1,i,\$2} } } }' \
        > .interval_list

        while read interval; do
            echo \$interval > \$(echo \${interval} | sed 's/[:*]/_/g' | sed 's/ /_/g').bed
        done < .interval_list

        """
}

process getGenomicIntervalList() {
    tag "extracting invervals from ${params.interval}..."
    cache 'lenient'
    output:
        path "*"
    script:
        """
        while read interval; do
            #check format of intervals
            if [[ -z \$(echo \${interval} | awk '{print \$3}') ]]; then
                echo \${interval} > \$(echo \${interval} | sed 's/[:*]/_/g').list
            else
                echo \${interval} > \$(echo \${interval} | sed 's/[:*]/_/g' | sed 's/ /_/g').bed
            fi
        done < ${params.interval}

        """
}

process haplotypeCaller() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'variantCaller'
    cache 'lenient'
    storeDir "${params.output_dir}/gvcfs/"
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        tuple \
            val(bamName), \
            path("${bamName}.g.vcf.gz"), \
            path("${bamName}.g.vcf.gz.tbi")
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            HaplotypeCaller \
            -I ${bamFile} \
            -R ${params.fastaRef} \
            -ERC GVCF \
            -OBI 'false' \
            -O ${bamName}.g.vcf.gz
        """
}

process haplotypeCallerWithIntervals() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'variantCaller'
    cache 'lenient'
    publishDir \
        path: "${params.output_dir}/gvcfs/intervals/"
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex), \
            path(interval)
    output:
        tuple \
            val(bamName), \
            path("${bamName}_${interval.simpleName}.g.vcf.gz"), \
            path("${bamName}_${interval.simpleName}.g.vcf.gz.tbi"), \
            path(interval)
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            HaplotypeCaller \
            -I ${bamFile} \
            -R ${params.fastaRef} \
            -L ${interval} \
            -ERC GVCF \
            -OBI 'false' \
            -O ${bamName}_${interval.simpleName}.g.vcf.gz
        """
}


// UNDER DEVELOPMENT
process concatGvcfs() {
    tag "processing ${bamName}"
    label 'bcftools'
    label 'variantCaller'
    cache 'lenient'
    //publishDir \
    //    path: "${params.output_dir}/gvcfs/", \
    //    mode: 'copy'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex), \
            path(interval)
    output:
        tuple \
            val(bamName), \
            path("${bamName}_${interval.simpleName}.g.vcf.gz"), \
            path("${bamName}_${interval.simpleName}.g.vcf.gz.tbi"), \
            path(interval)
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            HaplotypeCaller \
            -I ${bamFile} \
            -R ${params.fastaRef} \
            -L ${interval} \
            -ERC GVCF \
            -OBI 'false' \
            -O ${bamName}_${interval.simpleName}.g.vcf.gz
        """
}

process haplotypeCallerSpark() {
    tag "processing ${bamName}"
    beforeScript 'ulimit -c unlimited'
    label 'gatk'
    label 'variantCaller'
    cache 'lenient'
    publishDir \
        path: "${params.output_dir}/gvcfs/"
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        path "${bamName}.g.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            HaplotypeCallerSpark \
            -I ${bamFile} \
            -R ${params.fastaRef} \
            --native-pair-hmm-threads ${task.cpus} \
            -ERC GVCF \
            -OBI 'false' \
            -O "${bamName}.g.vcf.gz" \
            -- \
            --spark-master local[${task.cpus}]
        """
}

process combineGvcfs() {
    tag "Combining GVCFs to ${params.output_prefix}.gvcf.gz"
    label 'gatk' 
    label 'variantCaller'
    publishDir \
        path: "${params.output_dir}/gvcfs/combined/"
    input:
        path gvcfList
    output:
        path "${params.output_prefix}.g.vcf.{gz,gz.tbi}"
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "-V \${file}";
            fi
        done > gvcf.list

        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            CombineGVCFs \
            -R ${params.fastaRef} \
            --arguments_file gvcf.list \
            --dbsnp ${params.dbsnp} \
            -O "${params.output_prefix}.g.vcf.gz"
        """
}

process combineGvcfsPerChromosome() {
    tag "Combining GVCFs to ${params.output_prefix}.gvcf.gz"
    label 'gatk'
    label 'variantCaller'
    publishDir \
        path: "${params.output_dir}/gvcfs/combined/"
    input:
        tuple \
            val(interval), \
            path(gvcfList)
    output:
        tuple \
            val(interval), \
            path("${interval}_${params.output_prefix}.g.vcf.{gz,gz.tbi}")
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "-V \${file}";
            fi
        done > gvcf.list

        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            CombineGVCFs \
            -R ${params.fastaRef} \
            --arguments_file gvcf.list \
            -L ${interval} \
            --dbsnp ${params.dbsnp} \
            -O "${interval}_${params.output_prefix}.g.vcf.{gz,gz.tbi}"
        """
}

process createGenomicsDb() {
    tag "Combining GVCFs to ${params.output_prefix}.gvcf.gz"
    label 'gatk'
    label 'variantCaller'
    cache 'lenient'
    publishDir \
        path: "${params.output_dir}/genomicsdbs/"
    input:
        path gvcfList
    output:
        path "${params.output_prefix}-workspace"
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "-V \${file}";
            fi
        done > gvcf.list

        vcf=\$(ls *.g.vcf.gz | head -1)

        zgrep 'contig' \$(readlink \${vcf}) | \
            cut -f1 -d',' | \
            cut -f3 -d'=' > intervals.list

        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenomicsDBImport \
            -R ${params.fastaRef} \
            --tmp-dir . \
            --consolidate true \
            --arguments_file gvcf.list \
            -L intervals.list \
            --genomicsdb-workspace-path ${params.output_prefix}-workspace
        """
}

process createGenomicsDbPerInterval() {
    tag "processing ${interval.simpleName}..."
    label 'gatk'
    label 'genomisDBImport'
    publishDir \
        path: "${params.output_dir}/genomicsdbs/links/"      
    input:
        path(interval)
        path(gvcfList)
    output:
        tuple \
            val("${interval.simpleName}"), \
            path("${interval.simpleName}_${params.output_prefix}-workspace")
    script:
        """
        #for file in ${gvcfList}; do
        #    if [ \${file##*.} != "tbi" ]; then
        #        echo "-V \${file}";
        #    fi
        #done > gvcf.list
        mkdir -p temp
        gatk \
            --java-options "-Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ConcGCThreads=${task.cpus} -XX:ParallelGCThreads=${task.cpus}" \
            GenomicsDBImport \
            -R ${params.fastaRef} \
            --tmp-dir temp \
            --batch-size ${params.batch_size} \
            --consolidate true \
            --arguments_file ${gvcfList} \
            -L ${interval} \
            --genomicsdb-workspace-path ${interval.simpleName}_${params.output_prefix}-workspace
        """
}

//-XX:ConcGCThreads=${task.cpus} -XX:ParallelGCThreads=${task.cpus}

process updateGenomicsDbPerInterval() {
    tag "processing ${interval.simpleName}..."
    label 'gatk'
    label 'genomisDBImport'
    input:
        tuple \
            val(workspaceName), \
            path(interval), \
            path(workspace)
        path(gvcfList)
    output:
        tuple \
            val("${interval.simpleName}"), \
            path("${interval.simpleName}_${params.output_prefix}-workspace")
    script:
        """
        mkdir -p temp
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenomicsDBImport \
            -R ${params.fastaRef} \
            --tmp-dir temp \
            --consolidate true \
            --arguments_file ${gvcfList} \
            --batch-size ${params.batch_size} \
            -L ${interval} \
            --genomicsdb-update-workspace-path ${workspace}
        """
}

process getPedFile() {
    tag "Writing PED file"
    label 'bcftools'
    label 'smallMemory'
    echo true
    input:
        tuple \
            path(gvcf), \
            path(gvcfIndex)
    output:
        publishDir path: "${params.output_dir}/ped/"
        path "${params.output_prefix}.tmp.ped"
    script:
        """
        if [ ${params.ped} == NULL ]; then
           bcftools query -l ${gvcf} | \
               awk '{print \$1,\$1,"0","0","-9","-9"}' \
               > ${params.output_prefix}.ped
           ped="${params.output_prefix}.ped"
        else
           ped="${params.ped}"
        fi
        cp \${ped} ${params.output_prefix}.tmp.ped
        """
}

process genotypeGvcfs() {
    tag "Writing genotypes to ${params.output_prefix}.vcf.gz"
    label 'gatk'
    label 'variantCaller'
    publishDir \
        path: "${params.output_dir}/vcf/"
    input:
        tuple \
            path(gvcf), \
            path(gvcfIndex), \
            path(ped)
    output:
        path "${params.output_prefix}.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenotypeGVCFs \
            -R ${params.fastaRef} \
            --dbsnp ${params.dbsnp} \
            -ped ${ped} \
            -V ${gvcf} \
            -O "${params.output_prefix}.vcf.gz"
        """
}

process callVariantsFromGenomicsDB() {
    tag "Writing genotypes to ${interval.simpleName}_${params.output_prefix}.vcf.gz"
    label 'gatk'
    label 'variantCaller'
    //publishDir \
    //    path: "${params.output_dir}/vcf/"
    input:
        tuple \
            val(workspaceName), \
            path(interval), \
            path(workspace)
    output:
        path "${interval.simpleName}_${params.output_prefix}.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenotypeGVCFs \
            -R ${params.fastaRef} \
            --dbsnp ${params.dbsnp} \
            -L ${interval} \
            -V gendb://${workspace} \
            -O "${interval.simpleName}_${params.output_prefix}.vcf.gz"
        """
}

process callVariantsFromExistingGenomicsDB() {
    tag "Writing genotypes to ${workspaceName}.vcf.gz"
    label 'gatk'
    label 'variantCaller'
    //publishDir \
    //    path: "${params.output_dir}/vcf/"
    input:
        tuple \
            val(workspaceName), \
            path(workspace)
    output:
        path "${workspaceName}.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenotypeGVCFs \
            -R ${params.fastaRef} \
            --dbsnp ${params.dbsnp} \
            -V gendb://${workspace} \
            -O "${workspaceName}.vcf.gz"
        """
}

process collectIntervalsPerChromosome() {
    tag "Collecting intervals per chromosome..."
    label 'bcftools'
    label 'variantCaller'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    //storeDir "${params.output_dir}/vcf/"
    input:
        path(vcfs)
    output:
        path("*_vcfs_list.txt")
    script:
        """
        ls *.vcf.gz | awk '{print \$1,\$1}' > vcfs_list.txt

        while read line; do 
            data=( \$line );
            echo \$(basename \${data[0]} | sed 's/_/ /g' | awk '{print \$1,\$2}') \$(readlink \${data[1]})
        done < vcfs_list.txt | sort -g -k1 -k2 > vcf_chr_list.txt


        for chrom in \$(awk '{print \$1}' vcf_chr_list.txt | sort -V | uniq); do
            awk -v chr="\${chrom}" '\$1 == chr' vcf_chr_list.txt > \${chrom}_vcfs_list.txt
        done
        """
}

process concatPerChromIntervalVcfs() {
    tag "Concatenating VCF files per chromosome..."
    label 'bcftools'
    label 'variantCaller'
    //storeDir "${params.output_dir}/vcf/"
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    input:
        path(vcf_list)
    output:
        path("*_${params.output_prefix}.vcf.{gz,gz.tbi}")
    script:
        """
        chrom=\$(awk '{print \$1}' ${vcf_list} | uniq)
        awk '{print \$3}' ${vcf_list} > \${chrom}_concat.list

        bcftools \
            concat \
            -a \
            --threads ${task.cpus} \
            -Oz \
            -f \${chrom}_concat.list | \
        tee \${chrom}_${params.output_prefix}.vcf.gz | \
        bcftools \
            index \
            --threads ${task.cpus} \
            -ft \
            -o \${chrom}_${params.output_prefix}.vcf.gz.tbi
        """
}

process concatPerChromosomeVcfs() {
    tag "Concatenating all VCF files into ${params.output_prefix}.vcf.gz..."
    label 'bcftools'
    label 'variantCaller'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'move'
    input:
        path(vcf_list)
    output:
        path("${params.output_prefix}.vcf.{gz,gz.tbi}")
    script:
        """
        readlink *.gz | sort -V > concat.list

        bcftools \
            concat \
            -a \
            --threads ${task.cpus} \
            -Oz \
            -f concat.list | \
        tee ${params.output_prefix}.vcf.gz | \
        bcftools \
            index \
            --threads ${task.cpus} \
            -ft \
            -o ${params.output_prefix}.vcf.gz.tbi
        """
}

process deepVariantCaller() {
    tag "Writing genotypes to ${params.output_prefix}.vcf.gz"
    label 'deepvariant'
    label 'deepv_caller'
    beforeScript = 'module load python/anaconda-python-3.7'
    publishDir \
        path: "${params.output_dir}/gvcfs/"
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        path "${bamName}.g.vcf.{gz,gz.tbi}"
    script:
        """
        run_deepvariant \
            --model_type=\$([ ${params.exome} == true ] && echo WES || echo WGS) \
            --ref=${params.fastaRef} \
            --reads=${bamFile} \
            --output_gvcf=${bamName}.g.vcf.gz \
            --output_vcf=${bamName}.vcf.gz \
            --num_shards=${task.cpus}
        """
}

process glnexusJointCaller() {
    tag "Writing genotypes to ${params.output_prefix}.vcf.gz"
    label 'glnexus'
    label 'nexus_caller'
    publishDir \
        path: "${params.output_dir}/vcf/bcf/"
    input:
        path gvcfList
    output:
        path "${params.output_prefix}_glnexus.bcf"
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "\${file}";
            fi
        done | sort -V > gvcf.list

        if [ -d "GLnexus.DB" ]; then rm -rf GLnexus.DB; fi

        glnexus_cli \
            --config DeepVariant \
            --list gvcf.list \
            --threads ${task.cpus} \
            > "${params.output_prefix}_glnexus.bcf"
        """
}

process convertBcfToVcf() {
    tag "Writing genotypes to ${bcf_file.baseName}.vcf.gz"
    label 'bcftools'
    label 'nexus_caller'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    input:
        path bcf_file
    output:
        path "${bcf_file.baseName}.vcf.gz"
    script:
        """
        bcftools \
            view \
            --threads ${task.cpus} \
            -Oz \
            -o ${bcf_file.baseName}.vcf.gz \
            ${bcf_file}
        """
}

process dysguCallSvs() {
    tag "Writing genotypes to ${bamName}.vcf.gz"
    label 'dysgu'
    label 'dysgu_caller'
    publishDir \
        path: "${params.output_dir}/vcf/dysgu/", \
        mode: 'copy'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        path "${bamName}.vcf.gz"
    script:
        """
        dysgu \
            run \
            -p ${task.cpus} \
            ${params.fastaRef} \
            -x temp \
            --keep-small \
            --clean \
            -v2 \
            ${bamFile} | \
            bgzip -c > ${bamName}.vcf.gz \
        """
}

process indexVcf() {
    tag "processing ${vcf}"
    label 'bcftools'
    label 'mediumMemory'
    publishDir \
        path: "${params.output_dir}/vcf/dysgu/", \
        mode: 'copy'
    input:
        path(vcf)
    output:
        tuple \
            path(vcf), \
            path("${vcf}.tbi")
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${vcf}
        """
}

process dysguMergeVcfs() {
    tag "Writing genotypes to ${params.output_prefix}_dysgu_sv.vcf.gz"
    label 'dysgu'
    label 'dysgu_caller'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    input:
        path(vcfs)
    output:
        path "${params.output_prefix}_dysgu_sv.vcf.gz"
        """
        mkdir -p temp
        dysgu \
          merge \
          *.vcf.gz \
          --wd temp/ \
          --clean \
          -p ${task.cpus} | \
          bgzip -c > ${params.output_prefix}_dysgu_sv.vcf.gz
        """
}

process getBamChuncks() {
    tag "Splitting list of BAM files into chunks of 5..."
    input:
        path bam
    output:
        path("${params.output_prefix}_bamchunk_aaaaaaa*.txt")
    script:
        """
        readlink \$(ls *.{bam,cram}) > ${params.output_prefix}_bamlist.txt
        split \
            -l 5 \
            -a 8 \
            --additional-suffix .txt \
            ${params.output_prefix}_bamlist.txt \
            ${params.output_prefix}_bamchunk_
        """
}

process mantaCallSvs() {
    tag "Running MANTA..."
    label 'manta'
    label 'dysgu_caller'
    publishDir \
        path: "${params.output_dir}/vcf/manta/", \
        mode: 'copy'
    input:
        path bamlist
    output:
        path("*.vcf*")
    script:
        """
        bamlist=${bamlist}
        cat ${bamlist} | \
            awk '{print "--bam",\$1,"\\\\"}' | \
            sed "1 i --referenceFasta ${params.fastaRef} \$(if [[ ${params.exome} == true ]]; then echo --exome; fi)" | \
            sed 's/--exome/\\n--exome \\\\/g' | \
            sort -g | \
            sed '1 i configManta.py \\\\' \
            > \${bamlist/.txt/.sh}

        chmod 755 \${bamlist/.txt/.sh}

        ./\${bamlist/.txt/.sh}

        ./MantaWorkflow/runWorkflow.py \
            -j ${task.cpus}

        # ensure that there are no VCF files in the workdir
        [ -e "*.vcf*" ] && rm *.vcf*

        # get VCF results from manta workflow dir to nextflow workdir
        cp ./MantaWorkflow/results/variants/* .

        # rename manta results to include chunk filename
        for result in *.vcf*; do
            mv \${result} \${bamlist/.txt/}_\${result}
        done
        """
}

process mergeMantaDiploidSvCalls() {
    tag "mergeing diploid SV files..."
    label 'bcftools'
    label 'dysgu_caller'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    input:
        path manta_calls
    output:
        path("*.vcf.gz*")
    script:
        """
        bcftools \
            merge \
            -m all \
            --threads ${task.cpus} \
            -Oz \
            \$(ls *_diploidSV.vcf.gz | tr '\n' ' ') | \
            tee ${params.output_prefix}_manta_diploidSV.vcf.gz | \
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            -o ${params.output_prefix}_manta_diploidSV.vcf.gz.tbi
        """
}

process mergeMantaCandidateSmallIndelCalls() {
    tag "mergeing candidate small indels files..."
    label 'bcftools'
    label 'dysgu_caller'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    input:
        path manta_calls
    output:
        path("*.vcf.gz*")
    script:
        """
        bcftools \
            merge \
            -m all \
            --threads ${task.cpus} \
            -Oz \
            \$(ls *_candidateSmallIndels.vcf.gz | tr '\n' ' ') | \
            tee ${params.output_prefix}_manta_candidateSmallIndels.vcf.gz | \
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            -o ${params.output_prefix}_manta_candidateSmallIndels.vcf.gz.tbi
        """
}

process mergeMantaCandidateSvCalls() {
    tag "mergeing candidate SV files..."
    label 'bcftools'
    label 'dysgu_caller'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    input:
        path manta_calls
    output:
        path("*.vcf.gz*")
    script:
        """
        bcftools \
            merge \
            -m all \
            --threads ${task.cpus} \
            -Oz \
            \$(ls *_candidateSV.vcf.gz | tr '\n' ' ') | \
            tee ${params.output_prefix}_manta_CandidateSvs.vcf.gz | \
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            -o ${params.output_prefix}_manta_CandidateSvs.vcf.gz.tbi
        """
}
