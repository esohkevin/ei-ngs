#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//nextflow.enable.moduleBinaries = true

include {
    getInputFastqs;
    getSEInputFastqs;
    getInputAlignments;
    getAlignment;
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
        alignment = getAlignment().view()
    }

    if(!(params.buildVersion == 't2t')) {
        recalTable = recalibrateBaseQualityScores(alignment)
        alignment.combine(recalTable, by: 0).set { applyBQSR_input }
        recalibrated = applyBaseQualityRecalibrator(applyBQSR_input)
    }
    else {
        recalTable = recalibrateBaseQualityScores(alignment)
        alignment.combine(recalTable, by: 0).set { applyBQSR_input }
        recalibrated = applyBaseQualityRecalibrator(applyBQSR_input)
    }

}

workflow.onComplete { 
    println "\nDone! Check results in ${params.output_dir}\n" 
}
