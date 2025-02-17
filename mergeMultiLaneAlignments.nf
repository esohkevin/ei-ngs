#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//nextflow.enable.moduleBinaries = true

include {
    getInputFastqs;
    getSEInputFastqs;
    getAlignmentDir;
    sortAlignmentByName;
    sortCram;
    indexAlignment;
    indexAndMoveAlignment;
    fixAlignmentTags;
    fixAlignmentMate;
    mergeMultiLaneAlignment;
    updateMergedAlignmentHeader;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nMERGE ALIGNMENT FILES\n"
    getAlignmentDir()
        .set { alignment }
    mergeMultiLaneAlignment(alignment)
        .set { merged_alignment }
    updateMergedAlignmentHeader(merged_alignment)
        .set { merged_alignment_updated }
//    indexAndMoveAlignment(merged_alignment_updated)
}

workflow.onComplete { 
    println "\nDone! Check results in ${params.output_dir}\n" 
}
