#!/usr/bin/env python

#    do quality trimming for both reads with --nextseq-trim 2-color quality trimming
#    note: the minimum length requirement (trimmed read >= 5nt) is to sanitize the output for bowtie 1
rule clean_reads:
    input:
        r1 = lambda wc: SAMPLES[wc.sample]["r1"],
        r2 = lambda wc: SAMPLES[wc.sample]["r2"],
    output:
        r1 = "fastq/cleaned/{sample}-cleaned.r1.fastq.gz",
        r2 = "fastq/cleaned/{sample}-cleaned.r2.fastq.gz",
        log = "logs/clean_reads/clean_reads-{sample}.log"
    params:
        qual_cutoff = config["cutadapt"]["qual_cutoff"],
    conda:
        "../envs/cutadapt.yaml"
    threads:
        config["threads"]
    shell: """
        (cutadapt --cores={threads} --nextseq-trim={params.qual_cutoff} --minimum-length=5 --output={output.r1} --paired-output={output.r2} {input.r1} {input.r2}) &> {output.log}
        """

