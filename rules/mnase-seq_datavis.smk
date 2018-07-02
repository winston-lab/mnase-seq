#!/usr/bin/env python

localrules: cat_matrices

rule compute_matrix:
    input:
        annotation = lambda wc: FIGURES[wc.figure]["annotations"][wc.annotation]["path"],
        bw = "coverage/{norm}/{sample}_mnase-{readtype}-{norm}.bw"
    output:
        dtfile = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{readtype}-{norm}.mat.gz"),
        matrix = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{readtype}-{norm}.tsv"),
        melted = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{readtype}-{norm}-melted.tsv.gz")
    params:
        group = lambda wc : SAMPLES[wc.sample]["group"],
        refpoint = lambda wc: "TSS" if FIGURES[wc.figure]["parameters"]["type"]=="scaled" else FIGURES[wc.figure]["parameters"]["refpoint"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        binsize = lambda wc: FIGURES[wc.figure]["parameters"]["binsize"],
        binstat = lambda wc: FIGURES[wc.figure]["parameters"]["binstat"],
        nan_afterend = lambda wc: "--nanAfterEnd" if FIGURES[wc.figure]["parameters"]["nan_afterend"] and FIGURES[wc.figure]["parameters"]["type"]=="absolute" else [],
        anno_label = lambda wc: FIGURES[wc.figure]["annotations"][wc.annotation]["label"]
    threads : config["threads"]
    log: "logs/compute_matrix/compute_matrix-{figure}_{annotation}_{sample}_{readtype}_{norm}.log"
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute":
            shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} {params.nan_afterend} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {params.scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {params.refpoint} -g {params.group} -s {wildcards.sample} -a {params.anno_label} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wc: expand("datavis/{figure}/{norm}/{annotation}_{sample}-{readtype}-{norm}-melted.tsv.gz", annotation=[k for k,v in FIGURES[wc.figure]["annotations"].items()], sample=SAMPLES, figure=wc.figure, norm=wc.norm, readtype=wc.readtype)
    output:
        "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-{readtype}-{norm}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{figure}-{readtype}-{norm}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_figures:
    input:
        matrix = "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-{readtype}-{norm}.tsv.gz",
        annotations = lambda wc: [v["path"] for k,v in FIGURES[wc.figure]["annotations"].items()]
    output:
        heatmap_sample = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-seq_{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bysample.svg",
        heatmap_group = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-seq_{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bygroup.svg",
        meta_sample = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-seq_{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bysample.svg",
        meta_sample_overlay = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-seq_{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bysample-overlay.svg",
        meta_group = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-seq_{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bygroup.svg",
        meta_sample_clust = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-seq_{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bycluster-sample.svg",
        meta_group_clust = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-seq_{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bycluster-group.svg",
        metahmap_sample = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-seq_{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metahmap-bysample.svg",
        metahmap_group = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-seq_{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metahmap-bygroup.svg",
    params:
        # abusing snakemake a bit here...using params as output paths since in order to use lambda functions
        annotations_out = lambda wc: ["datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/".format(**wc) + annotation + "_cluster-" + str(cluster) + ".bed" for annotation in FIGURES[wc.figure]["annotations"] for cluster in range(1, FIGURES[wc.figure]["annotations"][annotation]["n_clusters"]+1)],
        clusters_out = lambda wc: ["datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/".format(**wc) + annotation + ".pdf" for annotation in FIGURES[wc.figure]["annotations"]],
        samplelist = lambda wc: get_samples(wc.status, wc.norm, [wc.condition, wc.control]),
        plottype = lambda wc: FIGURES[wc.figure]["parameters"]["type"],
        readtype = lambda wc: "dyad signal" if wc.readtype=="midpoint" else "protection",
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        pct_cutoff = lambda wc: FIGURES[wc.figure]["parameters"]["pct_cutoff"],
        spread_type = lambda wc: FIGURES[wc.figure]["parameters"]["spread_type"],
        trim_pct = lambda wc: FIGURES[wc.figure]["parameters"]["trim_pct"],
        refpointlabel = lambda wc: FIGURES[wc.figure]["parameters"]["refpointlabel"],
        endlabel = lambda wc:  "HAIL SATAN" if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["endlabel"],
        cmap = lambda wc: FIGURES[wc.figure]["parameters"]["heatmap_colormap"],
        sortmethod = lambda wc: FIGURES[wc.figure]["parameters"]["arrange"],
        cluster_scale = lambda wc: "FALSE" if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else str(FIGURES[wc.figure]["parameters"]["cluster_scale"]).upper(),
        cluster_samples = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else get_samples(wc.status, wc.norm, FIGURES[wc.figure]["parameters"]["cluster_conditions"]),
        cluster_five = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_five"],
        cluster_three = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_three"],
        k = lambda wc: [v["n_clusters"] for k,v in FIGURES[wc.figure]["annotations"].items()],
    conda: "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_mnase_figures.R"

