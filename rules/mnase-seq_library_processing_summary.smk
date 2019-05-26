#!/usr/bin/env python

localrules:
    build_read_processing_table,
    plot_read_processing,
    plot_spikein_pct,
    plot_fragment_lengths

rule get_fragment_lengths:
    input:
        expand("alignment/{sample}_mnase-seq-experimental.bam", sample=SAMPLES) if SISAMPLES else expand("alignment/{sample}_mnase-seq.bam", sample=SAMPLES)
    output:
        "qual_ctrl/fragment_length_distributions/mnase-seq_fragment_length_distributions.tsv"
    params:
        header = "\t".join(["fragsize"] + list(SAMPLES.keys()))
    threads: config["threads"]
    run:
        bam = input[0]
        shell("""samtools view {bam} | cut -f9 | sed 's/-//g' | sort -k1,1n -S 80% --parallel {threads} | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $2, $1}}' > {output}""")
        for bam in input[1:]:
            shell("""join -1 1 -2 2 -t $'\t' -e 0 -a 1 -a 2 --nocheck-order {output} <(samtools view {bam} | cut -f9 | sed 's/-//g' | sort -k1,1n -S 80% --parallel {threads} | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}') > qual_ctrl/fragment_length_distributions/.frag_length.temp; mv qual_ctrl/fragment_length_distributions/.frag_length.temp {output}""")
        shell("""sed -i "1i {params.header}" {output}""")

rule plot_fragment_lengths:
    input:
        table = "qual_ctrl/fragment_length_distributions/mnase-seq_fragment_length_distributions.tsv"
    output:
        plot = "qual_ctrl/fragment_length_distributions/mnase-seq_fragment_length_distributions.svg"
    conda: "../envs/tidyverse.yaml"
    script:
        "../scripts/mnase_frag_length.R"

rule build_read_processing_table:
    input:
        adapter = expand("logs/clean_reads/clean_reads-{sample}.log", sample=SAMPLES),
        align = expand("logs/align/align-{sample}.log", sample=SAMPLES),
    output:
        "qual_ctrl/read_processing/mnase-seq_read_processing_summary.tsv"
    log: "logs/read_processing_summary.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped" > {output}) &> {log}""")
        for sample, adapter, align in zip(SAMPLES.keys(), input.adapter, input.align):
            shell("""(grep -e "Total read pairs processed:" -e "Pairs written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &> {log}""")
            shell("""(grep -e "^Reported" {align} | awk '{{print $2}}' >> {output}) &> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing/mnase-seq_read_processing_summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing/mnase-seq_read_processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing/mnase-seq_read_processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing/mnase-seq_read_processing-loss.svg",
    conda: "../envs/tidyverse.yaml"
    script: "../scripts/processing_summary.R"

rule build_spikein_counts_table:
    input:
        total_bam = expand("alignment/{sample}_mnase-seq.bam", sample=SISAMPLES),
        exp_bam = expand("alignment/{sample}_mnase-seq-experimental.bam", sample=SISAMPLES),
        si_bam = expand("alignment/{sample}_mnase-seq-spikein.bam", sample=SISAMPLES),
    params:
        groups = [v["group"] for k,v in SISAMPLES.items()]
    output:
        "qual_ctrl/spikein/mnase-seq_spikein-counts.tsv"
    log:
        "logs/build_spikein_counts_table.log"
    run:
        shell("""(echo -e "sample\tgroup\ttotal_fragments\texperimental_fragments\tspikein_fragments" > {output}) &> {log} """)
        for sample, group, total_bam, exp_bam, si_bam in zip(SISAMPLES.keys(), params.groups, input.total_bam, input.exp_bam, input.si_bam):
            shell("""(paste <(echo -e "{sample}\t{group}") <(samtools view -c {total_bam}) <(samtools view -c {exp_bam}) <(samtools view -c {si_bam}) >> {output}) &>> {log}""")

rule plot_spikein_pct:
    input:
        "qual_ctrl/spikein/mnase-seq_spikein-counts.tsv"
    output:
        plot = "qual_ctrl/spikein/mnase-seq_spikein-plots-{status}.svg",
        stats = "qual_ctrl/spikein/mnase-seq_spikein-stats-{status}.tsv"
    params:
        samplelist = lambda wc : list(SISAMPLES.keys()) if wc.status=="all" else list(SIPASSING.keys()),
        conditions = conditiongroups_si if comparisons_si else [],
        controls = controlgroups_si if comparisons_si else [],
    conda: "../envs/tidyverse.yaml"
    script: "../scripts/plot_si_pct.R"

