#!/bin/bash

bsub -q priority -W 12:00 -n 1 -e snakemake.err -o snakemake.log snakemake -p --cluster-config cluster.yaml --use-conda --cluster "bsub -q {cluster.queue} -n {cluster.n} -W {cluster.time} -R 'rusage[mem={cluster.mem}]' -J {cluster.name} -eo {cluster.err} -oo {cluster.log}" --jobs 999
