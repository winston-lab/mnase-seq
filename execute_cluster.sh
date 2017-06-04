#!/bin/bash

bsub -q priority -W 24:00 -n 4 snakemake -p --cluster-config cluster.yaml --cluster "bsub -q {cluster.queue} -n {cluster.n} -W {cluster.time} -R 'rusage[mem={cluster.mem}]' -e {cluster.err} -o {cluster.out}" --jobs 50


