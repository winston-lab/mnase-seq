#!/bin/bash

bsub -q priority -W 24:00 -n 1 -e snakemake.err -o snakemake.log snakemake -p --use-conda --cluster-config cluster.yaml --cluster "bsub -q {cluster.queue} -n {cluster.n} -W {cluster.time} -R 'rusage[mem={cluster.mem}]' -J {cluster.name} -eo {cluster.err} -oo {cluster.log}" --jobs 25


