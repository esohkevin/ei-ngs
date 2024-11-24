#!/usr/bin/bash

#--- genemapgwas workflow wrapper ---#

projectdir=$(echo $(dirname ${0}))

ANSIRESET='\e[0m'
ANSIRED='\e[31m'
ANSIGRN='\e[32m'
ANSIYLW='\e[33m'
ANSIBLU='\e[34m'
ANSIPPL='\e[35m'
ANSIGRY='\e[2m'


function banner() {
   echo -e """
   ${ANSIGRY}===================================================================${ANSIRESET}
   ${ANSIBLU}Gene${ANSIRED}MAP${ANSIPPL}-NGS${ANSIRESET} ~ ${ANSIGRN}a wrapper for the nextflow-based genemapngs workflow${ANSIRESET}
   ${ANSIGRY}===================================================================${ANSIRESET}
   """
}

function checkprojectparams() {
   #projectdir=$(echo $(dirname ${0}))
   if [ ! -e ${projectdir}/nextflow.config ]; then
      echo -e "\n${ANSIRED}ERROR${ANSIRESET}: '${projectdir}/nextflow.config' not found!"
      echo -e "Please generate one from ${projectdir}/system.config ."
      echo -e "See workflow documentation on usage.\n"      
      exit 1
   else
      grep '=' ${projectdir}/nextflow.config | \
	 grep \
	   -e 'account' \
	   -e 'queue' \
	   -e 'containers_dir' \
	   -e 'workspace' \
         > ${projectdir}/.nextflow.config

      while read line; do
        param=( $(echo $line | grep '=' | sed "s|'||g") )
        paramname=${param[0]}
	paramval=$(echo ${param[2]}| tr [:lower:] [:upper:])
        if [[ ${paramval} == NULL ]] || [[ ${paramval} == "" ]]; then
          echo -e "\n${ANSIRED}ERROR${ANSIRESET}: ${paramname} = '${paramval}' ... FAIL!"
          echo -e "Please set a value for ${paramname} in '${projectdir}/nextflow.config'.\n"
          exit 1
        fi
      done < ${projectdir}/.nextflow.config
      rm ${projectdir}/.nextflow.config

      pn=( $(grep '=' ${projectdir}/nextflow.config | grep -w 'project_name' 2>/dev/null | sed "s|'||g") )
 
      if [[ ${pn[2]} == NULL ]] || [[ ${pn[2]} == "" ]]; then
         msg="""
      ${ANSIYLW}NOTE${ANSIRESET}:
      ${paramname} in '${ANSIGRN}${projectdir}/nextflow.config${ANSIRESET}' is not set! It has been initialized to
      '${ANSIGRN}myproject${ANSIRESET}'. This will be used as basename for all configuration files and some
      output files. If you wish to use a different name, edit '${ANSIGRN}${projectdir}/nextflow.config${ANSIRESET}'."""
      
         noteruler="""
      $( printf '+%.0s' {1..7} )"""
      
         echo -en "${noteruler}"
         echo -en "${msg}"
         echo -en "${noteruler}\n"
         echo ""
      
         sed -i -e "/project_name/s/.*/   project_name = 'myproject'/" ${projectdir}/nextflow.config
         projectname=myproject
         sleep 1
      else
         projectname=${pn[2]}
      
         msg="""
      ${ANSIYLW}NOTE${ANSIRESET}: 
      The project name '${ANSIGRN}${projectname}${ANSIRESET}' in '${ANSIGRN}${projectdir}/nextflow.config${ANSIRESET}' 
      will be used as basename for all configuration files and some output files.
      If you wish to use a different name, edit '${ANSIGRN}${projectdir}/nextflow.config${ANSIRESET}'."""
      
         noteruler="""
      $( printf '+%.0s' {1..7} )"""
      
         echo -en "${noteruler}"
         echo -en "${msg}"
         echo -en "${noteruler}\n"
         echo ""
         sleep 1
      fi
  fi
}

# display genemapngs banner
banner

# check and set project name in master config file
checkprojectparams

function usage() {
   echo -e """
   Usage: genemapngs <workflow> <profile> [options] ...

           workflows:
           ---------
               test: Run test to see if workfow installed correctly.
                 qc: Check FASTQ or Alignment (BAM/CRAM) quality.
               trim: Trim adapters and poor quality bases from reads.
              align: Align/map reads to reference and post-alignment processing.
         mergealign: Megre Alignment (BAM/CRAM) files.
            varcall: Perform variant calling (both single and joint sample) in one run.
           svarcall: Perform only sinlge sample variant calling to generate gVCF files.
           jvarcall: Perform only joint (multi-sample) variant calling with pre-existing gVCF files.
          varfilter: Filter variant calls in VCF/BCF files.


           profiles: <executor>,<container>,<reference>
           ---------
          executors: local, slurm
         containers: singularity, aptainer, docker
          reference: hg19, hg38, t2t


           examples:
           ---------
         genemapngs align slurm,singularity,hg38 [options]
         genemapngs qc local,singularity,hg19 --ftype FASTQ --output_dir ~/my-input-path/ --threads 1 --input_dir ~/my-output-path/ --njobs 1
   """
}

profiles() {
   echo """
           profiles: <executor>,<container>,<reference>
           ---------
          executors: local, slurm, pbs(pro)
         containers: singularity, aptainer, docker
          reference: hg19, hg38, t2t

	       e.g.: local,docker,hg19
	             slurm,singularity,hg38

   """
}

#####################################################################################################################
function checkprofile() {
   #profile was passed as only argument 
   #so it takes position $1 here
   if [[ "$1" == "" ]]; then
      echo "ERROR: please specify a profile to use!";
      profiles;
      exit 1;
   elif [[ $1 == -* ]]; then
      echo "ERROR: please specify a valid profile!";
      profiles;
      exit 1;
   else
      local profile="$1"
   fi
}

function check_ftype() {
   ftype=$1
   if [[ $ftype == -* ]]; then
      echo "ERROR: Invalid paramter value for option '--ftype'"
      exit 1;
   fi
}

function check_required_params() {
   for params_vals in $@; do
      #get each param and its value as an array
      param_val=( $(echo ${params_vals} | sed 's/,/ /g') )
      
      #slice the array to its consituent params and values
      param=${param_val[0]}
      val=${param_val[1]}

      #now check each param and its value
      if [[ $val == NULL ]]; then
         echo "ERROR: Missing required parameter '--${param}'";
         exit 1;
      elif [[ $val == -* ]]; then
         echo "ERROR: Invalid paramter value for option '--${param}'";
         exit 1;
      fi
   done
}

function check_optional_params() {
   for params_vals in $@; do
      #get each param and its value as an array
      param_val=( $(echo ${params_vals} | sed 's/,/ /g') )

      #slice the array to its consituent params and values
      param=${param_val[0]}
      val=${param_val[1]}

      #now check each param and its value
      if [[ $val == -* ]]; then
         echo "ERROR: Invalid paramter value for option '--${param}'";
         #break;
         exit 1;
      fi
   done
}

function setglobalparams() {
#resource selector file
rselector=$1
#- create the project nextflow config file
echo """includeConfig \"\${projectDir}/nextflow.config\"
includeConfig \"\${projectDir}/configs/profile-selector.config\"
includeConfig \"\${projectDir}/configs/resourceselector/${rselector}\"
"""
}


##################################################### USAGE #########################################################
function qcusage() {
   echo -e "\nUsage: genemapngs qc <profile> [options] ..."
   echo """
           options:
           --------
           --wgs                : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
                                  This is important for resource allocation.
           --ftype              : Input file type; FASTQ, BAM, CRAM [default: FASTQ].
           --input_dir          : (required) Path to FASTQ/BAM/CRAM files.
           --output_dir         : (optional) Results will be saved to parent of input directory ['input_dir/../'].
           --threads            : number of computer cpus to use [default: 4].
           --njobs              : (optional) number of jobs to submit at once [default: 4]
           --help               : print this help message.
   """
}

function trimusage() {
   echo -e "\nUsage: genemapngs trim <profile> [options] ..."
   echo """
           options:
           --------
           --wgs                : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
                                  This is important for resource allocation.
           --ftype              : Input file type; FASTQ, BAM, CRAM [default: FASTQ].
           --input_dir          : (required).
           --output_dir         : (optional) Results will be saved to parent of input directory ['input_dir/../'].
           --trimmer            : (optional) Options: trimmomatic, trimgalore [default: trimgalore].
           --adapter            : required if 'trimmomatic' selected [default: NP].
 
                         options:
                                  NP   --> NexteraPE-PE.fa
                                  T3U  --> TruSeq3-PE-2.fa [Illumina universal]
                                  T2P  --> TruSeq2-PE.fa
                                  T2S  --> TruSeq2-SE.fa
                                  T3P  --> TruSeq3-PE.fa
                                  T3S  --> TruSeq3-SE.fa

                              NB: 'trimgalore' will auto-detect adapters. Hence suitable to process
                                  samples from different sequencing companies in one batch.

           --min_length         : minimum read leangth to keep [default: 36].
           --headcrop           : number of bases to remove from the start of reads [default: 5].
           --crop               : number of bases to remove from the end of reads [default: 5].
           --threads            : number of computer cpus to use [default: 8].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}

function alignusage() {
   echo -e "\nUsage: genemapngs align <profile> [options] ..."
   echo """
           options:
           --------
           
           --wgs                : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
                                  This is important for resource allocation.
           --ftype              : Input file type; FASTQ, BAM, CRAM [default: FASTQ].
           --input_dir          : (required) Path to FASTQ/BAM/CRAM files.
           --aligner            : Alignment tool (required); BWA, DRAGMAP [default: BWA].
           --output_dir         : (optional) [results will be saved to parent of input directory].
           --se                 : If FASTQs are single-end (paired-end is assumed by default).
           --dup_marker         : Duplicate marker tool (optional); sambamba, samtools [default: sambamba]
           --remove_dup         : Whether to remove duplicates (optional): true, false [default: false]
           --spark              : Whether to use GATK spark mode for post-alignment processing (it multi-threads). Is not used by default.
           --threads            : number of computer cpus to use  [default: 11].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}

function mergealignusage() {
   echo -e "\nUsage: genemapngs align <profile> [options] ..."
   echo """
           options:
           --------
           
           --wgs                : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
                                  This is important for resource allocation.
           --input_dir          : (required) Path containing BAM/CRAM sub-directories.
                                  Each sub-directory contains only the files to be merged. The name of the sub-directory is used
                                  to name the resulting merged aligment file (CRAM).
           --output_dir         : (optional) [results will be saved to parent of input directory].
           --threads            : number of computer cpus to use  [default: 11].
           --sort_order         : (optional) options are: name, coordinate [default: coordinate]
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}

function varcallusage() {
   echo -e "\nUsage: genemapngs varcall <profile> [options] ..."
   echo """
           options:
           --------

           --wgs                : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
                                  This is important for resource allocation.
           --alignment_dir      : (required) Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).
           --out                : Output prefix (optional) [default: my-ngs].
           --output_dir         : (optional) [results will be saved to parent of input directory]
           --scaller            : Single sample variant caller; gatk, deepvariant [default: gatk]
                                  For structural variant calling, use 'svarcall'
           --jcaller            : Joint sample variant caller; gatk, glnexus [default: gatk]
           --batch_size         : (optional) number of samples to read into memory by GATK sample reader per time [default: 50]
           --interval           : List containing genomic intervals, one chromosome name per line and/or coordinate 
                                  in bed format: <chr> <start> <stop>.
                                  NB: Ensure that your chromosome names are the same as in the reference (e.g. chr1 or 1).
                                  If not provided, intervals will be created from gVCF header.
           --threads            : number of computer cpus to use  [default: 11].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}

function svarcallusage() {
   echo -e "\nUsage: genemapngs svarcall <profile> [options] ..."
   echo """
           options:
           --------

           --wgs                : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
                                  This is important for resource allocation.
           --alignment_dir      : (required) Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).
           --out                : Output prefix (optional) [default: my-ngs].
           --output_dir         : (optional) [results will be saved to parent of input directory]
           --scaller            : Single sample variant caller; gatk, deepvariant, dysgu, manta [default: gatk]
           --threads            : number of computer cpus to use  [default: 11].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}

function jvarcallusage() {
   echo -e "\nUsage: genemapngs jvarcall <profile> [options] ..."
   echo """
           options:
           --------

           --wgs                       : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
                                         This is important for resource allocation.
           --genomicsdb_workspace_dir  : (required) for calling variants from, or imprting gVCFS to, an existing genomicsdb workspace.
                                         Path containing genomicsdb workspaces (workspaces are directories).
           --imprt                     : Specify this flag if importing gVCF files to new genomicsdb workspaces. Cannot be used with '--update'
           --update                    : Specify this flag if importing gVCF files to existing genomicsdb workspaces. Cannot be used with '--imprt'
           --gvcf_dir                  : (required if '--update' flag is set). Path containing gVCF files and their indexes ('.tbi').
           --batch_size                : (optional) number of samples to read into memory by GATK sample reader per time [default: 50]
           --out                       : Output prefix (optional) [default: my-ngs].
           --output_dir                : (optional) [results will be saved to parent of input directory]
           --jcaller                   : Joint sample variant caller; gatk, glnexus [default: gatk]
           --interval                  : (required if '--imprt' or '--update' is specified) list containing genomic intervals.
                                         E.g. one chromosome name per line and/or coordinate in bed format: <chr> <start> <stop>.
                                         NB: Ensure that your chromosome names are the same as in the reference (e.g. chr1).
                                         If not provided, intervals will be created from gVCF header. If '--update' is specified and
                                         interval list is provided, it must be the same as that used for imprt.
           --threads                   : number of computer cpus to use  [default: 11].
           --njobs                     : (optional) number of jobs to submit at once [default: 10]
           --help                      : print this help message.
   """
}

function varfilterusage() {
   echo -e "\nUsage: genemapngs varfilter <profile> [options] ..."
   echo """
           options:
           --------

           --wgs                : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
                                  This is important for resource allocation.
	   --left_norm          : (optional) Whether to left-normalized variants such as is recommended by ANNOVAR; true, false [defaul: false].
           --vcf_dir            : (required) Path to VCF file(s).
           --minDP              : Minimum allele depth [default: 10].
           --minGQ              : Minimun genotype quality [default: 20].
           --minAC              : Minimun allele count (to remove singletons, set to 2) [default: 1].
           --out                : Output prefix (optional) [default: my-varfilter].
           --output_dir         : (optional) [results will be saved to parent of input directory]
           --jcaller            : The tool used to perform joint variant calling; gatk, glnexus [default: gatk]
           --threads            : number of computer cpus to use  [default: 4].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}


############################################# CONFIGURATION FILES ####################################################
function testconfig() {
#check and remove test config file if it exists
[ -e test.config ] && rm test.config

# $indir $bpm $csv $cluster $fasta $bam $out $outdir $thrds
echo """includeConfig \"\${projectDir}/nextflow.config\"
includeConfig \"\${projectDir}/configs/profile-selector.config\"
includeConfig \"\${projectDir}/configs/test.config\"
""" >> test.config
}


function qcconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${projectname}-qc.config ] && rm ${projectname}-qc.config

#qcconfig $exome $input_ftype $input_dir $output_dir $threads $njobs
echo """
params {
  //=======================================================
  //genemapngs reads quality check (qc) workflow parameters
  //=======================================================

  exome = $1
  input_ftype = '$2'
  input_dir = '$3'
  output_dir = '$4'
  threads = ${5}
  njobs = ${6}


  /*****************************************************************************************
  ~ exome: (optional) for GLNexus variant calling, manta structural variant calling, and for 
    resource management.
  ~ input_ftype: (required) FASTQ, BAM, CRAM.
  ~ input_dir: (required) Path to FASTQ/BAM/CRAM files.
  ~ output_dir: (optional) defaults to parent of input directory ['input_dir/../'].
  ~ threads: (optional) number of computer cpus to use  [default: 8]
  ~ njobs: (optional) number of jobs to submit at once [default: 10]
  *******************************************************************************************/
}

$(setglobalparams ${7})
""" >> ${projectname}-qc.config

echo -e "configuration file '${projectname}-qc.config' created!\n"
}


function trimconfig() { #params passed as arguments

#check and remove config file if it exists
[ -e ${projectname}-trim.config ] && rm ${projectname}-trim.config

#trimconfig $ftype $input_dir $output_dir $trimmer $adapter $min_length $headcrop $crop $threads $njobs
echo """
params {
  //====================================================
  //genemapngs reads trimming (trim) workflow parameters
  //====================================================

  exome = $1
  input_ftype = '$2'
  input_dir = '$3'
  output_dir = '$4'
  trimmer = '$5'
  adapter = '$6'
  min_length = ${7}
  headcrop = ${8}
  crop = ${9}
  threads = ${10}
  njobs = ${11}


  /*****************************************************************************************
  ~ exome: (optional) for GLNexus variant calling, manta structural variant calling, and for 
    resource management.
  ~ input_ftype: (required) FASTQ, BAM, CRAM.
  ~ input_dir: (required) Path to FASTQ/BAM/CRAM files.
  ~ output_dir: (optional) defaults to parent of input directory ['input_dir/../'].
  ~ trimmer: (optional) trimmomatic, trimgalore [default: trimgalore].
  ~ adapter: required if 'trimmomatic' selected [default: NP].
  ~ min_length: (optional) minimum read leangth to keep [default: 36].
  ~ headcrop: (optional) number of bases to remove from the start of reads [default: 5].
  ~ crop: (optional) number of bases to remove from the end of reads [default: 5].
  ~ threads: (optional) number of computer cpus to use  [default: 8]
  ~ njobs: (optional) number of jobs to submit at once [default: 10]
  *******************************************************************************************/
}

$(setglobalparams ${12})
""" >> ${projectname}-trim.config

echo -e "configuration file '${projectname}-trim.config' created!\n"
}


function alignconfig() { #params passed as arguments

#check and remove config file if it exists
[ -e ${projectname}-alignment.config ] && rm ${projectname}-alignment.config

#alignconfig $pe $exome $aligner $ftype $input_dir $output_dir $dup_marker $remove_dup $spark $threads $njobs
echo """
params {
  //=================================================
  // genemapngs alignment (align) workflow parameters 
  //=================================================

  pe = $1
  exome = $2
  aligner = '$3'
  input_ftype = '$4'
  input_dir = '$5'
  output_dir = '$6'
  dup_marker = '$7'
  remove_dup = $8
  spark = $9
  threads = ${10}
  njobs = ${11}


  /*****************************************************************************************
  ~ pe: (optional) whether reads are paired-end or single end): true, false [dfault: true].
  ~ exome: (optional) for GLNexus variant calling, manta structural variant calling, and for 
    resource management.
  ~ aligner: (optional) BWA, DRAGMAP [default: BWA].
  ~ input_ftype: (required) FASTQ, BAM, CRAM [default: FASTQ].
  ~ input_dir: (required) Path to FASTQ/BAM/CRAM files.
  ~ output_dir: (optional) defaults to parent of input directory ['input_dir/../']
  ~ dup_marker: (optional) duplicate marker tool: sambamba, samtools [default: sambamba]
    NB: This version (latest ~ https://hub.docker.com/r/maulik23/sambamba/tags) of sambamba 
    does not support CRAM (BAM outputs are generated). So samtools implemetation will conserve 
    space slightly better.
  ~ remove_dup: (optional) whether to remove duplicates true, false [default: false]
  ~ spark: (optional) use GATK spark mode for multi-threaded post-alignment processing
    true, false [default: false]
  ~ threads: (optional) number of computer cpus to use  [default: 11]
  ~ njobs: (optional) number of jobs to submit at once [default: 10]
  *******************************************************************************************/
}

$(setglobalparams ${12})
""" >> ${projectname}-alignment.config

echo -e "configuration file '${projectname}-alignment.config' created!\n"
}

function mergealignconfig() { #params passed as arguments

#check and remove config file if it exists
[ -e ${projectname}-merge-alignment.config ] && rm ${projectname}-merge-alignment.config

#mergealignconfig $exome $input_dir $output_dir $sort_order $threads $njobs
echo """
params {
  //======================================================
  // genemapngs alignment (mergealign) workflow parameters 
  //======================================================

  exome = $1
  input_dir = '$2'
  output_dir = '$3'
  sort_order = '${4}'
  threads = ${5}
  njobs = ${6}


  /*****************************************************************************************
  ~ exome: (optional) for GLNexus variant calling, manta structural variant calling, and for 
    resource management.
  ~ input_dir: (required) Path containing BAM/CRAM sub-directories.
    Each sub-directory contains only the files to be merged. The name of the sub-directory
    is used to name the resulting merged aligment file (CRAM).
  ~ output_dir: (optional) defaults to parent of input directory ['input_dir/../']
  ~ sort_order: (optional) options are: name, coordinate [default: coordinate]
  ~ threads: (optional) number of computer cpus to use  [default: 11]
  ~ njobs: (optional) number of jobs to submit at once [default: 10]
  *******************************************************************************************/
}

$(setglobalparams ${7})
""" >> ${projectname}-merge-alignment.config

echo -e "configuration file '${projectname}-merge-alignment.config' created!\n"
}

function varcallconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${projectname}-varcall.config ] && rm ${projectname}-varcall.config

#varcallconfig $exome $input_dir $output_dir $output_prefix $scaller $jcaller $threads $njobs
echo """
params {
  //=================================================================================
  // genemapngs single sample and joint variant calling (varcall) workflow parameters 
  //=================================================================================

  mode = 'varcall'
  exome = $1
  alignment_dir = '$2'
  output_dir = '$3'
  output_prefix = '$4'
  single_caller = '$5'
  joint_caller = '$6'
  batch_size = ${7}
  interval = '${8}'
  threads = ${9}
  njobs = ${10}

  
  /*****************************************************************************************
  ~ exome: (optional) for GLNexus variant calling, manta structural variant calling, and for 
    resource management.
  ~ alignment_dir: (required) path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).
  ~ output_dir: (optional) defaults to parent of input directory ['input_dir/../']
  ~ output_prefix: (optional) project name.
  ~ single_caller: (optional) gatk, deepvariant [default: gatk]
  ~ joint_caller: (optional) gatk, glnexus [default: gatk]
  ~ batch_size: (optional) number of samples to read into memory by GATK sample reader per
    time [default: 50]
  ~ interval: List containing genomic intervals, one chromosome name per line and/or coordinate 
    in bed format: <chr> <start> <stop>.
    NB: Ensure that your chromosome names are the same as in the reference (e.g. chr1 or 1).
    If not provided, intervals will be created from gVCF header.
  ~ threads: (optional) number of computer cpus to use  [default: 11]
  ~ njobs: (optional) number of jobs to submit at once [default: 10]
  *******************************************************************************************/   
}

$(setglobalparams ${11})
""" >> ${projectname}-varcall.config

echo -e "configuration file '${projectname}-varcall.config' created!\n"
}


function svarcallconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${projectname}-svarcall.config ] && rm ${projectname}-svarcall.config

#svarcallconfig $exome $alignment_dir $output_dir $output_prefix $scaller $threads $njobs
echo """
params {
  //=======================================================================
  //genemapngs single sample variant calling (svarcall) workflow parameters
  //=======================================================================

  mode = 'svarcall'
  exome = $1
  alignment_dir = '$2'
  output_dir = '$3'
  output_prefix = '$4'
  single_caller = '$5'
  threads = ${6}
  njobs = ${7}


  /*****************************************************************************************
  ~ exome: (optional) for GLNexus variant calling, manta structural variant calling, and for 
    resource management.
  ~ alignment_dir: (required) path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).
  ~ output_dir: (optional) defaults to parent of input directory ['input_dir/../']
  ~ output_prefix: (optional) project name.
  ~ single_caller: (optional) gatk, deepvariant, dysgu, manta [default: gatk]
  ~ threads: (optional) number of computer cpus to use  [default: 11]
  ~ njobs: (optional) number of jobs to submit at once [default: 10]
  *******************************************************************************************/ 
}

$(setglobalparams ${8})
""" >> ${projectname}-svarcall.config

echo -e "configuration file '${projectname}-svarcall.config' created!\n"
}


function jvarcallconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${projectname}-jvarcall.config ] && rm ${projectname}-jvarcall.config

#jvarcallconfig $exome $genomicsdb_workspace_dir $imprt $update $gvcf_dir $batch_size $output_dir $output_prefix $jcaller $interval $threads $njobs
echo """
params {
  //===============================================================
  //genemapngs joint variant calling (jvarcall) workflow parameters
  //===============================================================

  mode = 'jvarcall'
  exome = $1
  genomicsdb_workspace_dir = '$2'
  imprt = $3
  update = ${4}
  gvcf_dir = '$5'
  batch_size = ${6}
  output_dir = '$7'
  output_prefix = '$8'
  joint_caller = '$9'
  interval = '${10}'
  threads = ${11}
  njobs = ${12}


  /*****************************************************************************************
  ~ exome: (optional) for manta structural variant calling and resource management
  ~ genomicsdb_workspace_dir: (required) for calling variants from, or imprting gVCFS to, 
    an existing genomicsdb workspace. Path containing genomicsdb workspaces 
    (workspaces are directories).
  ~ imprt: (optional) whether to add gVCFs to NEW genomicsdb workspaces. 
    If true, 'gvcf_dir' must be provided [defaul: false]. Must be false if 'update = true'
  ~ update: (optional) whether to add gVCFs to EXISTING genomicsdb workspaces. 
    If true, 'gvcf_dir' must be provided [defaul: false]. Cannot be false if 'imprt = true'
  ~ gvcf_dir: required if 'update = true'. Path containing gVCF files and their indexes ('.tbi').
  ~ batch_size: (optional) number of samples to read into memory by GATK sample reader per 
    time [default: 50]
  ~ output_dir: (optional) defaults to parent of input directory 
    ['genomicsdb_workspace_dir/../' or 'gvcf_dir/../']
  ~ output_prefix: (required) name to add to output files.
  ~ joint_caller: (optional) gatk, deepvariant [default: gatk]
  ~ interval: (required if 'imprt' or 'update' is true) list containing genomic intervals. 
    E.g. one chromosome name per line and/or coordinate in bed format: <chr> <start> <stop>.
    NB: Ensure that your chromosome names are the same as in the reference (e.g. chr1).
    If not provided, intervals will be created from CRAM/gVCF header. If 'update = true' and 
    interval list is provided, it must be the same as that used for imprt.
  ~ threads: (optional) number of computer cpus to use  [default: 11]
  ~ njobs: (optional) number of jobs to submit at once [default: 10]
  *******************************************************************************************/
}

$(setglobalparams ${13})
""" >> ${projectname}-jvarcall.config

echo -e "configuration file '${projectname}-jvarcall.config' created!\n"
}


function varfilterconfig() {
#varfilterconfig $exome $vcf_dir $minDP $minGQ $minAC $out $output_dir $threads $njobs

#check and remove config file if it exists
[ -e ${projectname}-varfilter.config ] && rm ${projectname}-varfilter.config

echo """
params {
  //============================================================
  //genemapngs variant filtering (varfilter) workflow parameters
  //============================================================

  exome = $1
  left_norm = $2
  vcf_dir = '${3}'
  minDP = ${4}
  minGQ = ${5}
  minAC = ${6}
  output_prefix = '${7}'
  output_dir = '${8}'
  joint_caller = '${9}'
  threads = ${10}
  njobs = ${11}


  /*****************************************************************************************
  ~ exome: (optional) for VQSR, MQ annotation will be excluded for exome data
  ~ left_norm: (optional) Whether to left-normalized variants such as is recommended by 
    ANNOVAR; true, false [defaul: false].
  ~ vcf_dir: (required) path to VCF file(s).
  ~ minDP: (optional) minimum allele depth [default: 10]. 
  ~ minGQ: (optional) minimun genotype quality [default: 20]. 
  ~ minAC: (optional) Minimun allele count (to remove singletons, set to 2) [default: 1]. 
  ~ output_prefix: (optional) output prefix [default: my-varfilter]. 
  ~ output_dir: (optional) defaults to parent of vcf directory ['vcf_dir/../']
  ~ joint_caller: (optional) gatk, glnexus [default: gatk]
  ~ threads: (optional) number of computer cpus to use  [default: 11]
  ~ njobs: (optional) number of jobs to submit at once [default: 10]
  *******************************************************************************************/
}

$(setglobalparams ${12})
""" >> ${projectname}-varfilter.config

echo -e "configuration file '${projectname}-varfilter.config' created!\n"
}



if [ $# -lt 1 ]; then
   usage; exit 1;
else
   case $1 in
      test)
         profile='local,singularity,hg19'
         testconfig
      ;;
      qc)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            qcusage;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,ftype:,input_dir:,output_dir:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #- defaults
         exome=true
         resource=resource-selector-wes.config
         ftype=FASTQ     
         input_dir=NULL  
         output_dir=NULL 
         threads=4
         njobs=10       

         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; resource=resource-selector-wgs.config; shift;;
               --ftype) ftype="$2"; shift 2;;
               --input_dir) input_dir="$2"; shift 2;;
               --output_dir) output_dir="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; 2> /dev/null; shift 2;;
               --help) shift; qcusage; 2> /dev/null; exit 1;;
               --) shift; 2> /dev/null; break;;
               *) shift; qcusage; 2> /dev/null; exit 1;;
            esac
            continue; shift;
         done

         #########################################################
         # output_dir is required to avoid writing in the home   #
         # directory here the workflow should be preferrably run #
         #########################################################
         check_required_params \
	     input_dir,$input_dir \
             output_dir,$output_dir && \
         check_optional_params \
	     ftype,$ftype \
	     threads,$threads \
	     njobs,$njobs && \
         qcconfig \
             $exome \
	     $ftype \
	     $input_dir \
	     $output_dir \
	     $threads \
	     $njobs \
             $resource
      ;;
      trim)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            trimusage;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,ftype:,input_dir:,output_dir:,trimmer:,adapter:,min_length:,headcrop:,crop:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #- defaults
         exome=true
         resource=resource-selector-wes.config
         ftype=FASTQ       
         input_dir=NULL    
         output_dir=NULL    
         trimmer=trimgalore 
         adapter=NP         
         min_length=36      
         headcrop=5         
         crop=5             
         threads=8
         njobs=10           

         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; resource=resource-selector-wgs.config; shift;;
               --ftype) ftype="$2"; shift 2;;
               --input_dir) input_dir="$2"; shift 2;;
               --output_dir) output_dir="$2"; shift 2;;
               --trimmer) trimmer="$2"; shift 2;;
               --adapter) adapter="$2"; shift 2;;
               --min_length) min_length="$2"; shift 2;;
               --headcrop) headcrop="$2"; shift 2;;
               --crop) crop="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; trimusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; trimusage; 1>&2; exit 1;;
            esac
            continue; shift;
         done

         check_required_params \
	     input_dir,$input_dir \
             output_dir,$output_dir \
	     $([[ "${trimmer}" == "trimmomatic" ]] && echo "adapter,\$adapter") && \
         check_optional_params \
	     ftype,$ftype \
	     trimmer,$trimmer \
	     min_length,$min_length \
	     headcrop,$headcrop \
	     crop,$crop \
	     threads,$threads \
	     njobs,$njobs && \
         trimconfig \
             $exome \
	     $ftype \
	     $input_dir \
	     $output_dir \
	     $trimmer \
	     $adapter \
	     $min_length \
	     $headcrop \
	     $crop \
	     $threads \
	     $njobs \
             $resource

         #echo `nextflow -c ${out}-qc.config run qualitycontrol.nf -profile $profile`
      ;;
      align)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            alignusage; 1>&2;
            exit 1;
         fi

         prog=`getopt -a --long "help,se,wgs,spark,aligner:,ftype:,input_dir:,output_dir:,dup_marker:,remove_dup:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #defaults
         pe=true             
         aligner=BWA         
         ftype=FASTQ         
         input_dir=NULL      
         output_dir=NULL     
         exome=true
         resource=resource-selector-wes.config
         dup_marker=sambamba 
         remove_dup=false    
         spark=false
         threads=11
         njobs=10            
         
         eval set -- "$prog"

         while true; do
            case $1 in
               --se) pe=false; shift;;
               --wgs) exome=false; resource=resource-selector-wgs.config; shift;;
               --spark) spark=true; shift;;
               --aligner) aligner="${2}"; shift 2;;
               --ftype) ftype="$2"; shift 2;;
               --input_dir) input_dir="$2"; shift 2;;
               --output_dir) output_dir="${2}"; shift 2;;
               --dup_marker) dup_marker="$2"; shift 2;;
               --remove_dup) remove_dup="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; alignusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; alignusage; 1>&2; exit 1;;
            esac
         done

         check_required_params \
	     input_dir,$input_dir \
             output_dir,$output_dir && \
         check_optional_params \
	     ftype,$ftype \
	     dup_marker,$dup_marker \
	     remove_dup,$remove_dup \
	     threads,$threads \
	     njobs,$njobs && \
         alignconfig \
	     $pe \
	     $exome \
	     $aligner \
	     $ftype \
	     $input_dir \
	     $output_dir \
	     $dup_marker \
	     $remove_dup \
	     $spark \
	     $threads \
	     $njobs \
             $resource
      ;;
      mergealign)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            mergealignusage; 1>&2;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,sort_order:,input_dir:,output_dir:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #defaults
         input_dir=NULL      
         output_dir=NULL     
         exome=true
         resource=resource-selector-wes.config
         sort_order=coordinate
         threads=11
         njobs=10            
         
         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; resource=resource-selector-wgs.config; shift;;
               --sort_order) sort_order=$2; shift 2;;
               --input_dir) input_dir="$2"; shift 2;;
               --output_dir) output_dir="${2}"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; mergealignusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; mergealignusage; 1>&2; exit 1;;
            esac
         done

         check_required_params \
	     input_dir,$input_dir \
             output_dir,$output_dir && \
         check_optional_params \
	     sort_order,$sort_order \
	     threads,$threads \
	     njobs,$njobs && \
         mergealignconfig \
	     $exome \
	     $input_dir \
	     $output_dir \
	     $sort_order \
	     $threads \
	     $njobs \
             $resource
      ;;
      varcall)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            varcallusage; 1>&2;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,batch_size:,interval:,alignment_dir:,output_dir:,out:,scaller:,jcaller:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #defaults
         #ftype=FASTQ        
         alignment_dir=NULL   
         output_dir=NULL
         batch_size=50
         interval=NULL
         output_prefix="myngs"
         ped=NULL             
         scaller=gatk         
         exome=true
         resource=resource-selector-wes.config
         jcaller=gatk         
         threads=11
         njobs=10             

         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; resource=resource-selector-wgs.config; shift;;
               #--ftype) ftype="$2"; shift 2;;
               --alignment_dir) alignment_dir="$2"; shift 2;;
               --batch_size) batch_size=$2; shift 2;;
               --interval) interval=$2; shift 2;;
               --output_dir) output_dir="${2}"; shift 2;;
               --out) output_prefix="$2"; shift 2;;
               --scaller) scaller="$2"; shift 2;;
               --jcaller) jcaller="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; varcallusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; varcallusage; 1>&2; exit 1;;
            esac
         done

         check_required_params \
             output_dir,$output_dir \
	     alignment_dir,$alignment_dir && \
         check_optional_params \
             batch_size,$batch_size \
             interval,$interval \
	     output_prefix,$output_prefix \
	     scaller,$scaller \
	     jcaller,$jcaller \
	     threads,$threads \
	     njobs,$njobs && \
         varcallconfig \
	     $exome \
	     $alignment_dir \
	     $output_dir \
	     $output_prefix \
	     $scaller \
	     $jcaller \
             $batch_size \
             $interval \
	     $threads \
	     $njobs \
             $resource
           

      ;;
      svarcall)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            svarcallusage; 1>&2;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,alignment_dir:,output_dir:,out:,scaller:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #defaults
         #ftype=FASTQ        
         alignment_dir=NULL   
         output_dir=NULL      
         output_prefix="myngs"
         ped=NULL             
         scaller=gatk         
         exome=true
         resource=resource-selector-wes.config
         threads=11
         njobs=10             

         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; resource=resource-selector-wgs.config; shift;;
               #--ftype) ftype="$2"; shift 2;;
               --alignment_dir) alignment_dir="$2"; shift 2;;
               --output_dir) output_dir="${2}"; shift 2;;
               --out) output_prefix="$2"; shift 2;;
               --scaller) scaller="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; svarcallusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; svarcallusage; 1>&2; exit 1;;
            esac
         done

         check_required_params \
             output_dir,$output_dir \
	     alignment_dir,$alignment_dir && \
         check_optional_params \
	     output_prefix,$output_prefix \
	     scaller,$scaller \
	     threads,$threads \
	     njobs,$njobs && \
         svarcallconfig \
	     $exome \
	     $alignment_dir \
	     $output_dir \
	     $output_prefix \
	     $scaller \
	     $threads \
	     $njobs \
             $resource

         #echo `nextflow -c ${out}-idat2vcf.config run idat2vcf.nf -profile $profile`

      ;;
      jvarcall)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            jvarcallusage; 1>&2;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,gvcf_dir:,imprt,update,genomicsdb_workspace_dir:,batch_size:,output_dir:,out:,jcaller:,interval:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #defaults
         #ftype=FASTQ                 
         output_dir=NULL              
         output_prefix="myngs"       
         ped=NULL                     
         exome=true
         resource=resource-selector-wes.config
         jcaller=gatk                 
         interval=NULL                
         gvcf_dir=NULL                
         update=false
         imprt=false  #import appears to be a nextflow reserved variable name
         genomicsdb_workspace_dir=NULL 
         batch_size=50                 
         threads=11
         njobs=10                      

         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; resource=resource-selector-wgs.config; shift;;
               --update) update=true; shift;;
               --imprt) imprt=true; shift;;
               --output_dir) output_dir="${2}"; shift 2;;
               --out) output_prefix="$2"; shift 2;;
               --jcaller) jcaller="$2"; shift 2;;
               --interval) interval="$2"; shift 2;;
               --gvcf_dir) gvcf_dir="$2"; shift 2;;
               --genomicsdb_workspace_dir) genomicsdb_workspace_dir="$2"; shift 2;;
               --batch_size) batch_size="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; jvarcallusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; jvarcallusage; 1>&2; exit 1;;
            esac
         done

         if [[ "${imprt}" == "true" ]]; then

            ##############################################
            # If importing to new genomicsdb workspace,  #
            # then only gvcf_dir is required and update  #
            # must be false                              #
            ##############################################

            update=false
            check_required_params \
                output_dir,$output_dir \
                gvcf_dir,$gvcf_dir && \
            check_optional_params \
                output_prefix,$output_prefix \
                jcaller,$jcaller \
                interval,$interval \
                threads,$threads \
                njobs,$njobs \
                batch_size,$batch_size && \
            jvarcallconfig \
                $exome \
                $genomicsdb_workspace_dir \
                $imprt \
                $update \
                $gvcf_dir \
                $batch_size \
                $output_dir \
                $output_prefix \
                $jcaller \
                $interval \
                $threads \
                $njobs \
                $resource

            #echo `nextflow -c ${out}-idat2vcf.config run idat2vcf.nf -profile $profile`

         elif [[ "${update}" == "true" ]]; then

            ##############################################
            # If updating existing genomicsdb workspace, #
            # then genomicsdb_workspace_dir and interval #
            # list are also required. Must be the same   #
            # list as the one used for importing. imprt  #
            # option must necessarily be false           #
            ##############################################

            imprt=false
            check_required_params \
                output_dir,$output_dir \
	        gvcf_dir,$gvcf_dir \
                interval,$interval \
		genomicsdb_workspace_dir,$genomicsdb_workspace_dir && \
            check_optional_params \
                output_prefix,$output_prefix \
                jcaller,$jcaller \
                wgs,$exome \
                threads,$threads \
                njobs,$njobs \
                batch_size,$batch_size && \
            jvarcallconfig \
                $exome \
                $genomicsdb_workspace_dir \
                $imprt \
                $update \
                $gvcf_dir \
                $batch_size \
                $output_dir \
                $output_prefix \
                $jcaller \
                $interval \
                $threads \
                $njobs \
                $resource

            #echo `nextflow -c ${out}-idat2vcf.config run idat2vcf.nf -profile $profile`

         else
            check_required_params \
                output_dir,$output_dir \
	        genomicsdb_workspace_dir,$genomicsdb_workspace_dir && \
            check_optional_params \
                gvcf_dir,$gvcf_dir \
                output_prefix,$output_prefix \
                jcaller,$jcaller \
                interval,$interval \
                wgs,$exome \
                threads,$threads \
                njobs,$njobs \
                batch_size,$batch_size && \
            jvarcallconfig \
                $exome \
                $genomicsdb_workspace_dir \
                $imprt \
                $update \
                $gvcf_dir \
                $batch_size \
                $output_dir \
                $output_prefix \
                $jcaller \
                $interval \
                $threads \
                $njobs \
                $resource

            #echo `nextflow -c ${out}-idat2vcf.config run idat2vcf.nf -profile $profile`

         fi
      ;;
      varfilter)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            varfilterusage;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,left_norm,vcf_dir:,minDP:,minGQ:,minAC:,out:,output_dir:,jcaller:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #- defaults
         exome=true
         resource=resource-selector-wes.config
	 left_norm=false
         vcf_dir=NULL       
         minDP=10           
         minGQ=20           
         minAC=1            
         out="my-varfilter" 
         output_dir=NULL    
         jcaller=gatk       
         threads=4          
         njobs=10           

         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; resource=resource-selector-wgs.config; shift;;
               --left_norm) exome=true; shift;;
               --vcf_dir) vcf_dir="$2"; shift 2;;
               --minDP) minDP="$2"; shift 2;;
               --minGQ) minGQ="$2"; shift 2;;
               --minAC) minAC="$2"; shift 2;;
               --out) out="$2"; shift 2;;
               --output_dir) output_dir="$2"; shift 2;;
               --jcaller) jcaller="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; varfilterusage; exit 1;;
               --) shift; break;;
               *) shift; varfilterusage; exit 1;;
            esac
            continue; shift;
         done

         check_required_params \
             output_dir,$output_dir \
             vcf_dir,$vcf_dir && \
         check_optional_params \
             out,$out \
             minDP,$minDP \
             minGQ,$minGQ \
             minAC,$minAC \
             jcaller,$jcaller \
             threads,$threads \
             njobs,$njobs && \
         varfilterconfig \
             $exome \
	     $left_norm \
	     $vcf_dir \
	     $minDP \
	     $minGQ \
	     $minAC \
	     $out \
	     $output_dir \
	     $jcaller \
	     $threads \
	     $njobs \
             $resource

         #echo `nextflow -c ${out}-qc.config run qualitycontrol.nf -profile $profile`

      ;;
      help) shift; usage; exit 1;;
      *) shift; usage; exit 1;;
   esac
fi


