#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {

  println "\nGeneMAP NGS WORKFLOW: TEST\n"

  println "pe=${params.pe}"
  println "aligner=${params.aligner}"
  println "ftype=${params.input_ftype}"
  println "input_dir=${params.input_dir}"
  println "output_dir=${params.output_dir}"
  println "output_prefix=${params.output_prefix}"
  println "single_caller=${params.single_caller}"
  println "exome=${params.exome}"
  println "joint_caller=${params.joint_caller}"
  println "gvcf_dir=${params.gvcf_dir}"
  println "spark=${params.spark}"
  println "threads=${params.threads}"
  println "njobs=${params.njobs}"
  println ""

  run_plink()
	
}



workflow.onComplete { 
  println "Workflow completed at: ${workflow.complete}"
  println "     Execution status: ${ workflow.success ? 'OK' : 'failed'}"
}

workflow.onError{
  println "workflow execution stopped with the following message: ${workflow.errorMessage}"
}

process run_plink() {

  // directives
  tag "processing ..."
  label 'plink'
  label 'test'
  debug true
  //echo true
  
  script:
    """
    echo -e "\nPLINK2 is used for the test as it is light-weight and easily pulled from docker hub\n"

    plink2 \
      --help \
      --file
    """
}
