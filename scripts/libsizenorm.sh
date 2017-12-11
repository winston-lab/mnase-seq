#!/bin/bash

## normalize two.bedgraph by million counts present in one.bedgraph,
## scaling coverage by a scalefactor
##
## usage: libsizenorm.sh one.bedgraph two.bedgraph scalefactor > two-normalized-by-one.bedgraph
##

awk -v sf=$3 'BEGIN{FS=OFS="\t";sum=0}{if(NR==FNR){sum+=$4}else{print $1, $2, $3, $4*1000000*sf/sum}}' $1 $2
