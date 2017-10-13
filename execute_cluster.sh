#!/bin/bash

if [[ ! -z $SLURM_CONF ]]
then
    #O2/SLURM
    sbatch -p priority -t 12:00:00 --mem 4G -n 1 -e snakemake.err -o snakemake.log --wrap='snakemake -p --latency-wait 60 --cluster-config cluster.yaml --use-conda --cluster "sbatch -p {cluster.queue} -n {cluster.n} -t {cluster.time} --mem {cluster.mem} -J {cluster.name} -e {cluster.err} -o {cluster.log}" --jobs 999'
else
    #orchestra/LSF
    bsub -q priority -W 12:00 -n 1 -e snakemake.err -o snakemake.log snakemake -p --cluster-config cluster.yaml --use-conda --cluster "bsub -q {cluster.queue} -n {cluster.n} -W {cluster.time} -R 'rusage[mem={cluster.mem}]' -J {cluster.name} -eo {cluster.err} -oo {cluster.log}" --jobs 999
fi
