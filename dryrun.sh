#!/bin/bash

snakemake -np -R `cat <(snakemake --lc --rerun-incomplete) <(snakemake --li --rerun-incomplete) <(snakemake --lp --rerun-incomplete) | sort -u` --rerun-incomplete --use-conda
