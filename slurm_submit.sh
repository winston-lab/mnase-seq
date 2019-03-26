#!/bin/bash

#SBATCH -p short
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=1400M
#SBATCH -c 1
#SBATCH -e snakemake.err
#SBATCH -o snakemake.log
#SBATCH -J MNase-seq-snakemake
#SBATCH --mail-type=ALL                    # ALL email notification type
#SBATCH --mail-user=cweiner@g.harvard.edu  # Email to which notifications will be sent
snakemake -p \
    -R `cat <(snakemake --lc --rerun-incomplete) \
            <(snakemake --li --rerun-incomplete) \
            <(snakemake --lp --rerun-incomplete) | sort -u` \
    --latency-wait 300 \
    --rerun-incomplete \
    --cluster-config $(grep -h annotation_workflow config.yaml | \
                        head -n 1 | \
                        cut -f1 --complement -d ":" | \
                        awk '{print $1}' | \
                        paste -d '' - <(echo cluster.yaml)) \
    --cluster-config cluster.yaml \
    --use-conda \
    --jobs 9999 \
    --restart-times 1 \
    --cluster "sbatch -p {cluster.queue} -c {cluster.n} -t {cluster.time} --mem-per-cpu={cluster.mem} -J {cluster.name} -e {cluster.err} -o {cluster.log} --parsable" \
    --cluster-status "bash slurm_status.sh"

