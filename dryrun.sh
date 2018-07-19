#!/bin/bash

snakemake -np -R `cat <(snakemake --lc) <(snakemake --li) <(snakemake --lp)` --rerun-incomplete --use-conda
