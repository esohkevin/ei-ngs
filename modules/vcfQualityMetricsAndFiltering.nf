def getVcf() {
    return channel.fromPath( params.vcf_dir + "*.vcf.gz" )
}

def getThousandGenomesReference() {
    return channel.fromFilePairs( params.ref_dir + "chr*1kg.phase3.v5a.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error "\nAn error occurred! Please check that the reference file and its index '.tbi' exist...\n" }
	   	  .map { chr, ref_file -> 
			 tuple( chr.replaceFirst(/chr/,""), ref_file.first(), ref_file.last())
	    	  }
}

process getVcfIndex() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        path input_vcf
    output:
        tuple \
            path("${input_vcf}"), \
            path("${input_vcf}.tbi")
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${input_vcf}
        """
}

process splitVcfs() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'vqsr'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        tuple \
            path("${input_vcf.simpleName}.snp.vcf.gz"), \
            path("${input_vcf.simpleName}.indel.vcf.gz")
    script:
        """
        gatk \
           --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
           SplitVcfs \
           -I ${input_vcf} \
           --STRICT false \
           --SNP_OUTPUT ${input_vcf.simpleName}.snp.vcf.gz \
           --INDEL_OUTPUT ${input_vcf.simpleName}.indel.vcf.gz
        """
}

process vqsrSnp() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'vqsr'
    publishDir \
        path: "${params.output_dir}/vqsr-tables/", \
        mode: 'copy'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        tuple \
            path(input_vcf), \
            path(vcf_index), \
            path("${input_vcf.simpleName}.recal"), \
            path("${input_vcf.simpleName}.recal.idx"), \
            path("${input_vcf.simpleName}.tranches")
    script:
        """
        gatk \
           --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
           VariantRecalibrator \
           -R ${params.fastaRef} \
           -tranche 100.0 -tranche 99.95 -tranche 99.9 \
           -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
           -tranche 95.0 -tranche 94.0 \
           -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
           -V ${input_vcf} \
           --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
           --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.snpRef} \
           --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omniRef} \
           --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
           -an QD \
           \$([ ${params.exome} == false ] && echo "-an DP") \
           -an MQ \
           -an MQRankSum \
           -an ReadPosRankSum \
           -an FS \
           -an SOR \
           -mode SNP \
           -O ${input_vcf.simpleName}.recal \
           --tranches-file ${input_vcf.simpleName}.tranches \
           --rscript-file ${input_vcf.simpleName}.plots.R
        """
}

process vqsrIndel() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'vqsr'
    publishDir \
        path: "${params.output_dir}/vqsr-tables/", \
        mode: 'copy'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        tuple \
            path(input_vcf), \
            path(vcf_index), \
            path("${input_vcf.simpleName}.recal"), \
            path("${input_vcf.simpleName}.recal.idx"), \
            path("${input_vcf.simpleName}.tranches")
    script:
        """
        gatk \
           --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
           VariantRecalibrator \
           -R ${params.fastaRef} \
           -V ${input_vcf} \
           -tranche 100.0 -tranche 99.95 -tranche 99.9 \
           -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
           -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
           -tranche 92.0 -tranche 91.0 -tranche 90.0 \
           --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
           --resource:1000G,known=false,training=true,truth=false,prior=12.0 ${params.indelsRef} \
           --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omniRef} \
           --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
           -an QD \
           \$([ ${params.exome} == false ] && echo "-an DP") \
           -an MQRankSum \
           -an ReadPosRankSum \
           -an FS \
           -an SOR \
           -mode INDEL \
           --max-gaussians 4 \
           -O ${input_vcf.simpleName}.recal \
           --tranches-file ${input_vcf.simpleName}.tranches \
           --rscript-file ${input_vcf.simpleName}.plots.R
        """
}

process applyVqsrSnp() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'applyVqsr'
    publishDir \
        path: "${params.output_dir}/vqsr/", \
        mode: 'copy'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index), \
            path(recal_table), \
            path(recal_index), \
            path(tranches)   
    output:
        tuple \
            val("${input_vcf.simpleName}"), \
            path("${input_vcf.simpleName}.snp.vqsr.vcf.gz"), \
            path("${input_vcf.simpleName}.snp.vqsr.vcf.gz.tbi")
    script:
        """
        gatk \
            --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            ApplyVQSR \
            -V ${input_vcf} \
            --recal-file ${recal_table} \
            --mode SNP \
            --truth-sensitivity-filter-level 99.9 \
            --create-output-variant-index true \
            --tranches-file ${tranches} \
            -R ${params.fastaRef} \
            -O ${input_vcf.simpleName}.snp.vqsr.vcf.gz
        """
}

process applyVqsrIndel() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'applyVqsr'
    publishDir \
        path: "${params.output_dir}/vqsr/", \
        mode: 'copy'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index), \
            path(recal_table), \
            path(recal_index), \
            path(tranches)   
    output:
        tuple \
            val("${input_vcf.simpleName}"), \
            path("${input_vcf.simpleName}.indel.vqsr.vcf.gz"), \
            path("${input_vcf.simpleName}.indel.vqsr.vcf.gz.tbi")
    script:
        """
        gatk \
            --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            ApplyVQSR \
            -V ${input_vcf} \
            --recal-file ${recal_table} \
            --mode INDEL \
            --truth-sensitivity-filter-level 99.9 \
            --create-output-variant-index true \
            --tranches-file ${tranches} \
            -R ${params.fastaRef} \
            -O ${input_vcf.simpleName}.indel.vqsr.vcf.gz
        """
}

process mergeVCFs() {
    tag "merging VCFs"
    label 'bcftools'
    label 'longRun'
    publishDir \
        path: "${params.output_dir}/vqsr/", \
        mode: 'copy'
    input:
        tuple \
            val(vcfbase), \
            path(snp_vcf), \
            path(snp_vcf_index), \
            path(indel_vcf), \
            path(indel_vcf_index)
    output:
        tuple \
            val(vcfbase), \
            path("${vcfbase}.snp.indel.vqsr.vcf.gz"), \
            path("${vcfbase}.snp.indel.vqsr.vcf.gz.tbi")
    script:
        """
        bcftools \
            concat \
            -a \
            -d all \
            --threads ${task.cpus} \
            -Oz \
            $snp_vcf \
            $indel_vcf | \
            tee ${vcfbase}.snp.indel.vqsr.vcf.gz | \
        bcftools index \
            --threads ${task.cpus} \
            -ft \
            --output "${vcfbase}.snp.indel.vqsr.vcf.gz.tbi"
        """
}

process filterGatkCalls() {
    tag "VCF supplied: ${vcfbase}"
    label 'bcftools'
    label 'longRun'
    publishDir \
        path: "${params.output_dir}/filtered/", \
        mode: 'copy'
    input:
        tuple \
            val(vcfbase), \
            path(vcf), \
            path(index)
    output:
        tuple \
            val("${vcfbase}"), \
            path("${vcfbase}.filtered.vcf.gz")
    script:
        """
        bcftools \
            view \
            -i \'FILTER=="PASS"\' \
            --threads ${task.cpus} \
            -Oz \
            ${vcf} | \
        tee "${vcfbase}.filtered.vcf.gz" | \
        bcftools \
            index \
            --threads ${task.cpus} \
            -ft \
            --output "${vcfbase}.filtered.vcf.gz.tbi"
        """
}


process filterGlnexusCalls() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    publishDir \
        path: "${params.output_dir}/filtered/", \
        mode: 'copy'
    input:
        path(input_vcf)
        path(vcf_index)
    output:
        tuple \
            path("${input_vcf.simpleName}.filtered.vcf.gz"), \
            path("${input_vcf.simpleName}.filtered.vcf.gz.tbi")
    script:
        """
        bcftools \
            view \
            -i \'FILTER!="MONOALLELIC"\' \
            --threads ${task.cpus} \
            ${input_vcf} | \
            bcftools \
                view \
                -i \'DP>=${params.minDP}\' \
                --threads ${task.cpus} | \
                bcftools \
                    view \
                    -i \'GQ>=${params.minGQ}\' \
                    --threads ${task.cpus} \
                    -Oz | \
                    tee "${input_vcf.simpleName}.filtered.vcf.gz" | \
                bcftools index \
                    --threads ${task.cpus} \
                    -ft \
                    --output "${input_vcf.simpleName}.filtered.vcf.gz.tbi"
        """
}

process splitMultiallelicSnvs() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        tuple \
            path("${input_vcf.simpleName}-tmp.vcf.gz"), \
            path("${input_vcf.simpleName}-tmp.vcf.gz.tbi")
    script:
        """
        bcftools \
            norm \
            -m -any \
            --threads ${task.cpus} \
            -Oz \
            ${input_vcf} | \
            tee "${input_vcf.simpleName}-tmp.vcf.gz" | \
        bcftools index \
            --threads ${task.cpus} \
            -ft \
            --output "${input_vcf.simpleName}-tmp.vcf.gz.tbi"            
        """
}

process leftnormalizeSnvs() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    publishDir \
        path: "${params.output_dir}/filtered/", \
        mode: 'copy'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        tuple \
            path("${input_vcf.simpleName}-filtered-leftnorm.vcf.gz"), \
            path("${input_vcf.simpleName}-filtered-leftnorm.vcf.gz.tbi")
    script:
        """
        bcftools \
            norm \
            -f ${params.fastaRef} \
            --threads ${task.cpus} \
            -Oz \
            ${input_vcf} | \
        bcftools \
            view \
            -c \$([[ ${params.minAC} == null ]] && echo 1 || echo ${params.minAC}) \
            --threads ${task.cpus} \
            -Oz | \
            tee "${input_vcf.simpleName}-filtered-leftnorm.vcf.gz" | \
        bcftools index \
            --threads ${task.cpus} \
            -ft \
            --output "${input_vcf.simpleName}-filtered-leftnorm.vcf.gz.tbi"
        """
}

process getCleanVcf() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    publishDir \
        path: "${params.output_dir}/filtered/", \
        mode: 'copy'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        path "${input_vcf.simpleName}-filtered-leftnorm-clean.vcf.{gz,gz.tbi}"
    script:
        """
        bcftools \
            view \
            --threads ${task.cpus} \
            -r \$(echo chr{1..22}, chrX | sed 's/[[:space:]]//g') \
            -Oz \
            "${input_vcf}" | 
            tee "${input_vcf.simpleName}-filtered-leftnorm-clean.vcf.gz" | \
            bcftools index \
            --threads ${task.cpus} \
            -ft \
            --output "${input_vcf.simpleName}-filtered-leftnorm-clean.vcf.gz.tbi"
        """
}

process getVcfStats() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    publishDir \
        path: "${params.output_dir}/vcfstats/", \
        mode: 'copy'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        path "${input_vcf.simpleName}.vcfstats.txt"
    script:
        """
        bcftools \
            stats \
            -F "${params.fastaRef}" \
            -s - \
            "${input_vcf}" > \
            "${input_vcf.simpleName}.vcfstats.txt"
        """
}

process plotVcfStats() {
    tag "Processing ${vcfstat}"
    label 'bcftools'
    label 'longRun'
    publishDir \
        path: "${params.output_dir}_${vcfstat.simpleName}", \
        mode: 'copy'
    input:
        path vcfstat
    output:
        path "./*"
    script:
        """
        mkdir -p "${params.output_dir}${vcfstat.simpleName}"
        plot-vcfstats \
           -p . \
           "${vcfstat}"
        """
}

