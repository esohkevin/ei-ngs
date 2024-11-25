#!/usr/bin/env bash

######################################################
# make non-overlapping intervals of N base pairs for #
# chromosomes larger than N base pairs. Otherwise,   #
# print out chromosome length as only interval.      #
# Alternate and HLA contigs can be excluded so only  #
# major contigs (chromosomes) are processed.         #
# For HLA typing, use the SNP2HLA or equivalent tool #
######################################################

if [ $# -lt 4 ]; then
  echo -e "Usage: get_interval_list.sh [gvcf-file] [interval-size] [no-alt <yes/no>] [output-prefix]"
  echo ""
  echo "no-alt: exclude all alternate contigs"
else
  gvcf=$1  #"/scratch/eshkev001/projects/wes/higenes/careni/hg38/gvcfs/OS.P001_92401_S25_val.bqsr.g.vcf.gz"
  intval=$2
  noalt=$(echo $3 | tr [:upper:] [:lower:])
  out=$4
  container=/scratch/eshkev001/containers/sickleinafrica-bcftools-1.11.img

  zgrep \
    '##contig' \
    ${gvcf} \
    > gvcf.hdr

  if [[ "${noalt}" == "yes" ]]; then
    sed 's/[=,>]/\t/g' gvcf.hdr | \
        cut -f3,5 | \
        awk '
          {
            if($1 ~ /^chr/ && length($1) <= 5) {
              print $1,$2
            } else if(!($1 ~ /^chr/) && length($1) <= 2) {
              print $1,$2
            }
          }' | \
        awk \
          -v intvl=${intval} '
            {
              if( $2 <= intvl ) {
                print $1,"0",$2
              } else{
                  for( i=0; i<=$2; i+=intvl ) {
                    if( i+(intvl-1) < $2 ) {
                      print $1,i,i+(intvl-1)
                    } else{
                        print $1,i,$2
                    }
                  }
              }
            }' \
    > ${out}_noalt_interval_list.txt
  else
    sed 's/[=,>]/\t/g' gvcf.hdr | \
        cut -f3,5 | \
        awk \
          -v intvl=${intval} '
            { 
              if( $2 <= intvl ) { 
                print $1,"0",$2 
              } else{ 
                  for( i=0; i<=$2; i+=intvl ) { 
                    if( i+(intvl-1) < $2 ) { 
          	    print $1,i,i+(intvl-1)
                    } else{ 
          	      print $1,i,$2 
                    } 
                  } 
              } 
            }' \
    > ${out}_withalt_interval_list.txt
  fi
  rm gvcf.hdr
fi

#while read interval; do
#    echo $interval > $(echo ${interval} | sed 's/[:*]/_/g' | sed 's/ /_/g').bed
#done < .interval_list

