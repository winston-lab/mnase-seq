#!/usr/bin/env python

localrules:
    map_counts_to_transcripts,
    combine_counts_over_transcripts

rule map_counts_to_transcripts:
    input:
        bed = lambda wc: config["genome"]["transcript_annotation"] if wc.species=="experimental" else config["spike_in"]["transcript_annotation"],
        bg = lambda wc: f"coverage/counts/{wc.sample}_mnase-midpoint-counts.bedgraph" if wc.species=="experimental" else f"coverage/sicounts/{wc.sample}_mnase-midpoint-sicounts.bedgraph"
    output:
        temp("diff_levels/{condition}-v-{control}/{sample}_{species}-counts-over-transcripts.tsv")
    log: "logs/map_counts_to_transcripts/map_counts_to_transcripts_{condition}-v-{control}_{sample}-{species}.log"
    shell: """
        (LC_COLLATE=C sort -k1,1 -k2,2n {input.bed} | bedtools map -a stdin -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{($6=="+") ? strand="plus" : strand="minus"; print $4"~"$1"-"strand"~"$2"~"$3, $7}}' &> {output}) &> {log}
        """

rule combine_counts_over_transcripts:
    input:
        lambda wc : ["diff_levels/{condition}-v-{control}/".format(**wc) + x + f"_{wc.species}-counts-over-transcripts.tsv" for x in get_samples("passing", "libsizenorm", [wc.control, wc.condition])]
    output:
        "diff_levels/{condition}-v-{control}/{condition}-v-{control}_{species}-allsample-counts-over-transcripts.tsv"
    params:
        n = lambda wc: 2*len(get_samples("passing", "libsizenorm", [wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples("passing", "libsizenorm", [wc.control, wc.condition]))
    log: "logs/combine_counts_over_transcripts/combine_counts_over_transcripts-{condition}-v-{control}-{species}.log"
    shell: """
        (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
        """

rule call_nuclevel_changes:
    input:
        expcounts = "diff_levels/{condition}-v-{control}/{condition}-v-{control}_experimental-allsample-counts-over-transcripts.tsv",
        sicounts = lambda wc: [] if wc.norm=="libsizenorm" else "diff_levels/{condition}-v-{control}/{condition}-v-{control}_spikein-allsample-counts-over-transcripts.tsv".format(**wc)
    output:
        results_all = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-results-{norm}-all.tsv",
        results_up = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-results-{norm}-up.tsv",
        results_down = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-results-{norm}-down.tsv",
        results_unch = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-results-{norm}-unchanged.tsv",
        bed_all = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-results-{norm}-all.bed",
        bed_up = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-results-{norm}-up.bed",
        bed_down = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-results-{norm}-down.bed",
        bed_unch = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-results-{norm}-unchanged.bed",
        normcounts = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-counts-sfnorm-{norm}.tsv",
        rldcounts = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-counts-rlog-{norm}.tsv",
        qcplots = "diff_levels/{condition}-v-{control}/{norm}/{condition}-v-{control}-mnase-seq-qcplots-{norm}.svg"
    params:
        samples = lambda wc : get_samples("passing", wc.norm, [wc.control, wc.condition]),
        groups = lambda wc : [PASSING[x]["group"] for x in get_samples("passing", wc.norm, [wc.control, wc.condition])],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"])
    conda: "../envs/diff_exp.yaml"
    script:
        "../scripts/call_de_transcripts.R"

