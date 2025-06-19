#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//nextflow.enable.moduleBinaries = true

include {
    getInputFastqs;
    getSEInputFastqs;
    getInputAlignments;
    sortAlignmentByName;
    convertAlignmentToFastq;
    bwaAligner;
    alignReadsBWA;
    dragenAligner;
    tmapAligner;
    convertSamToBam;
    convertBamToCram;
    sortAlignment;
    sortAlignmentToBam;
    sortCram;
    indexAlignment;
    indexAlignment as indexCram;
    indexAndCopyAlignment;
    markDuplicatesGatk;
    markDuplicates;
    markDupSambam;
    indexAlignment as indexMarkedAlignment;
    recalibrateBaseQualityScores;
    applyBaseQualityRecalibrator;
    markDuplicatesSpark;
    fixAlignmentTags;
    fixAlignmentMate;
    recalibrateBaseQualityScoresSpark;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nAlignment workflow begins here\n"
    if( params.input_ftype.toUpperCase() == "FASTQ" ) {
        println "INPUT FILE TYPE IS FASTQ\n"
        if(params.pe == false) {
            println "SINGLE END READS\n"
            fastq = getSEInputFastqs().view()
        }
        else {
            println "PAIRED END READS\n"
            fastq = getInputFastqs()
        }
    }
    else {
        println "INPUT FILE TYPE IS ALIGNMENT (BAM/CRAM)\n"
        alignment = getInputAlignments()
        alignmentSortedByName = sortAlignmentByName(alignment)
        fastq = convertAlignmentToFastq(alignmentSortedByName)
    }

    if( params.aligner == "DRAGMAP" ) {
        sam = dragenAligner(fastq)
    }
    else if( params.aligner == "TMAP" ) {
        sam = tmapAligner(fastq)
    }
    else {
        //sam = bwaAligner(fastq)
        dupsMarked = alignReadsBWA(fastq)
    }

    // TEST GATK BUNDLES FOR T2T WITH BQSR //

    dupsMarkedIndexed = indexMarkedAlignment(dupsMarked)
    recalTable = recalibrateBaseQualityScores(dupsMarkedIndexed)
    dupsMarkedIndexed.combine(recalTable, by: 0).set { applyBQSR_input }
    recalibrated = applyBaseQualityRecalibrator(applyBQSR_input)


//    if(!(params.buildVersion == 't2t')) {
//        dupsMarkedIndexed = indexMarkedAlignment(dupsMarked)
//        recalTable = recalibrateBaseQualityScores(dupsMarkedIndexed)
//        dupsMarkedIndexed.combine(recalTable, by: 0).set { applyBQSR_input }
//        recalibrated = applyBaseQualityRecalibrator(applyBQSR_input)
//    }
//    else {
//        dupsMarkedIndexed = indexAndCopyAlignment(dupsMarked)
//    }

/*
*    alignment = fixAlignmentMate(sam)
*
*    if(params.buildVersion == 't2t') {
*        if(params.spark == true) {
*            dupsMarked = markDuplicatesSpark(alignment)
*            fixedAlignment = fixAlignmentTags(dupsMarked)
*        }
*        else {
*
*            //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
*            // BQSR is not implemented for t2t reference //
*            // because there are no gatk bundles for yet //
*            //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
*
*            if(params.dup_marker.toUpperCase() == "SAMTOOLS") {
*              sortedAlignment = sortAlignment(alignment)
*              dupsMarked = markDuplicates(sortedAlignment)
*              dupsMarkedIndexed = indexAndCopyAlignment(dupsMarked)
*            }
*            else {
*              sortedAlignment = sortAlignmentToBam(alignment)
*              dupsMarked = markDupSambam(sortedAlignment)
*              dupsMarkedIndexed = indexAndCopyAlignment(dupsMarked)
*            }
*        }
*    } 
*    else {
*        if(params.spark == true) {
*            dupsMarked = markDuplicatesSpark(alignment)
*            fixedAlignment = fixAlignmentTags(dupsMarked)
*            recalibrated = recalibrateBaseQualityScoresSpark(fixedAlignment)
*            cram = convertBamToCram(recalibrated)
*        }
*        else {
*            if(params.dup_marker.toUpperCase() == "SAMTOOLS") {
*              sortedAlignment = sortAlignment(alignment)
*              dupsMarked = markDuplicates(sortedAlignment)
*            }
*            else {
*              sortedAlignment = sortAlignmentToBam(alignment)
*              dupsMarked = markDupSambam(sortedAlignment)
*            }
*
*            dupsMarkedIndexed = indexMarkedAlignment(dupsMarked)
*            recalTable = recalibrateBaseQualityScores(dupsMarkedIndexed)
*            dupsMarkedIndexed.combine(recalTable, by: 0).set { applyBQSR_input }
*            recalibrated = applyBaseQualityRecalibrator(applyBQSR_input)
*            //cram = convertBamToCram(recalibrated)
*        }
*    }
*
*/

}

workflow.onComplete { 
    println "\nDone! Check results in ${params.output_dir}\n" 
}
