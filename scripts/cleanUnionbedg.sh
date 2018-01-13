#!/bin/bash

awk 'BEGIN{FS=OFS="\t"} NR==1{ORS="\t"; print "name"; for(k=4;k<NF;k++) print $k; ORS="\n"; print $NF} {ORS="\t"; sum=0; for(i=4;i<=NF;i++) sum+=$i} sum>0{print $1"-"$2"-"$3; for(j=4;j<NF;j++) print $j; ORS="\n"; print $NF}'
