#!/usr/bin/env python

configfile: "config.yaml"

SAMPLES = config["samples"]
sisamples = {k:v for (k,v) in SAMPLES.items() if v["spikein"]=="y"}
PASSING = {k:v for (k,v) in SAMPLES.items() if v["pass-qc"] == "pass"}
sipassing = {k:v for (k,v) in PASSING.items() if v["spikein"] == "y"}

controlgroups = config["comparisons"]["libsizenorm"]["controls"]
conditiongroups = config["comparisons"]["libsizenorm"]["conditions"]

if sipassing:
    controlgroups_si = config["comparisons"]["spikenorm"]["controls"]
    conditiongroups_si = config["comparisons"]["spikenorm"]["conditions"]

COUNTTYPES = ["counts", "sicounts"] if sisamples else ["counts"]
NORMS = ["libsizenorm", "spikenorm"] if sisamples else ["libsizenorm"]

localrules: all,
            make_barcode_file,
            bowtie_build,
            samtools_index,
            # build_nucwave_input,
            # get_fragsizes,
            # cat_fragsizes,
            # make_fragments_table,
            # make_window_file,
            # cat_windows,
            # cat_perbase_depth,
            # group_bam_for_danpos,

rule all:
    input:
        #fastqc
        "qual_ctrl/fastqc/raw",
        expand("qual_ctrl/fastqc/{sample}/{sample}-clean.{read}_fastqc.html", sample=SAMPLES, read=["r1","r2"]),
        #demultiplex
        expand("fastq/{sample}.{read}.fastq.gz", sample=SAMPLES, read=["r1","r2"]),
        #alignment
        expand("alignment/{sample}.bam", sample=SAMPLES),
        expand("alignment/unaligned-{sample}_{r}.fastq.gz", sample=SAMPLES, r=[1,2]),
        #coverage
        expand("coverage/{counttype}/{sample}-mnase-midpoint-{counttype}.bedgraph", sample=SAMPLES, counttype=COUNTTYPES),
        expand("coverage/{norm}/{sample}-mnase-{readtype}-{norm}.bedgraph", norm=NORMS, sample=SAMPLES, readtype=["midpoint"]),
        #datavis
        expand(expand("datavis/{{annotation}}/spikenorm/mnase-{{annotation}}-spikenorm-{{status}}_{condition}-v-{control}-{{readtype}}-{{plot}}-bysample.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), annotation=config["annotations"], readtype=["midpoint"], status=["all","passing"], plot=["heatmap","metagene"]) + expand(expand("datavis/{{annotation}}/libsizenorm/mnase-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}-{{readtype}}-{{plot}}-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], readtype=["midpoint"], status=["all","passing"], plot=["heatmap", "metagene"]) if sisamples else expand(expand("datavis/{{annotation}}/libsizenorm/mnase-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}-{{readtype}}-{{plot}}-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], readtype=["midpoint"], status=["all","passing"], plot=["heatmap","metagene"])

def plotcorrsamples(wildcards):
    dd = SAMPLES if wildcards.status=="all" else PASSING
    if wildcards.condition=="all":
        if wildcards.norm=="libsizenorm": #condition==all,norm==lib
            return list(dd.keys())
        else: #condition==all,norm==spike
            return [k for k,v in dd.items() if v["spikein"]=="y"]
    elif wildcards.norm=="libsizenorm": #condition!=all;norm==lib
        return [k for k,v in dd.items() if v["group"]==wildcards.control or v["group"]==wildcards.condition]
    else: #condition!=all;norm==spike
        return [k for k,v in dd.items() if (v["group"]==wildcards.control or v["group"]==wildcards.condition) and v["spikein"]=="y"]

rule make_barcode_file:
    output:
        "fastq/barcodes.tsv"
    run:
        with open(output[0], "w") as out:
            for k,v in SAMPLES.items():
                out.write(k + '\t' + v["barcode"] + '\n')

#fastQC on raw sequencing data
rule fastqc_raw:
    input:
       r1 = config["fastq"]["r1"],
       r2 = config["fastq"]["r2"]
    output:
       "qual_ctrl/fastqc/raw"
    threads : config["threads"]
    log : "logs/fastqc/fastqc-raw.log"
    shell: """
        (mkdir -p {output}) &> {log}
        (fastqc -o {output} --noextract -t {threads} {input.r1} {input.r2}) &>> {log}
        """

#cutadapt doesn't demultiplex paired end, so we use a different script
#NOTE: as of v1.15, cutadapt now demultiplexes paired end, could try it out
# this would allow us to throw out reads where adapter isn't found in both reads
rule demultiplex:
    input:
        r1 = config["fastq"]["r1"],
        r2 = config["fastq"]["r2"],
        barcodes = "fastq/barcodes.tsv"
    output:
        r1 = expand("fastq/{sample}.r1.fastq.gz", sample=SAMPLES),
        r2 = expand("fastq/{sample}.r2.fastq.gz", sample=SAMPLES)
    log:
        "logs/demultiplex.log"
    shell: """
       (fastq-multx -B {input.barcodes} -b {input.r1} {input.r2} -o fastq/%.r1.fastq.gz -o fastq/%.r2.fastq.gz ) &> {log}
        """
# note: fastq-multx only remove the barcode on read 1. The barcode on read 2 is removed with cutadapt.

# cutadapt:
#    remove barcode from read 2
#    trim the 'A' tail off
#    do quality trimming for both reads
#        - ideally, we would use --nextseq-trim for 2-color quality trimming instead of -q
#            - however, --nextseq-trim currently doesn't trim read 2
#    note: the minimum length requirement (trimmed read >= 5nt) is to sanitize the output for bowtie 1
#    note: the maximum length requirement is discard reads in which the barcode isn't found in read 2
rule cutadapt:
    input:
        r1 = "fastq/{sample}.r1.fastq.gz",
        r2 = "fastq/{sample}.r2.fastq.gz"
    output:
        r1 = "fastq/cleaned/{sample}-clean.r1.fastq.gz",
        r2 = "fastq/cleaned/{sample}-clean.r2.fastq.gz"
    params:
        qual_cutoff = config["cutadapt"]["qual_cutoff"],
        adapter = lambda wildcards : SAMPLES[wildcards.sample]["barcode"]+"T",
        max_len = lambda wildcards: config["read-length"] - len(SAMPLES[wildcards.sample]["barcode"]+"T")
    # threads: config["threads"]
    log:
        "logs/cutadapt/cutadapt-{sample}.log"
    shell: """
        (cutadapt -e 0.15 -u 1 -G ^{params.adapter} -q {params.qual_cutoff} --minimum-length 5 --maximum-length {params.max_len} -o {output.r1} -p {output.r2} {input.r1} {input.r2}) &> {log}
        """

#fastQC on demultiplexed and cleaned reads
rule fastqc_cleaned:
    input:
        r1 = "fastq/cleaned/{sample}-clean.r1.fastq.gz",
        r2 = "fastq/cleaned/{sample}-clean.r2.fastq.gz"
    output:
        "qual_ctrl/fastqc/{sample}/{sample}-clean.r1_fastqc.html",
        "qual_ctrl/fastqc/{sample}/{sample}-clean.r2_fastqc.html"
    threads : config["threads"]
    log : "logs/fastqc/fastqc-cleaned-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.sample}) &> {log}
        (fastqc -o qual_ctrl/fastqc/{wildcards.sample} --noextract -t {threads} {input.r1} {input.r2}) &>> {log}
        """

rule bowtie_build:
    input:
        fasta = config["combinedgenome"]["fasta"] if sisamples else config["genome"]["fasta"]
    output:
        expand(config["bowtie"]["index-path"] + "/{{basename}}.{num}.ebwt", num=[1,2,3,4]),
        expand(config["bowtie"]["index-path"] + "/{{basename}}.rev.{num}.ebwt", num=[1,2]),
    params:
        idx_path = config["bowtie"]["index-path"],
        prefix = config["combinedgenome"]["experimental_prefix"]
    log:
        "logs/bowtie-build.log"
    run:
        if sisamples:
            shell("(bowtie-build {input.fasta} {params.idx_path}/{wildcards.basename}) &> {log}")
        else:
            shell("(sed -e 's/>/>{params.prefix}/g' {input.fasta} > .{params.prefix}.fa; bowtie-build .{params.prefix}.fa {params.idx_path}/{wildcards.basename}; rm .{params.prefix}.fa) &> {log}")

#align with Bowtie 1
#in Christine's paper, Burak uses -m 10 --best
rule bowtie:
    input:
        expand(config["bowtie"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".{num}.ebwt", num=[1,2,3,4]) if sisamples else expand(config["bowtie"]["index-path"] + "/" + config["genome"]["name"] + ".{num}.ebwt", num=[1,2,3,4]),
        expand(config["bowtie"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".rev.{num}.ebwt", num=[1,2]) if sisamples else expand(config["bowtie"]["index-path"] + "/" + config["genome"]["name"] + ".rev.{num}.ebwt", num=[1,2]),
        r1 = "fastq/cleaned/{sample}-clean.r1.fastq.gz",
        r2 = "fastq/cleaned/{sample}-clean.r2.fastq.gz"
    params:
        idx_path = config["bowtie"]["index-path"],
        basename = config["combinedgenome"]["name"] if sisamples else config["genome"]["name"],
        max_mismatch = config["bowtie"]["max_mismatch"],
        min_ins = config["bowtie"]["min_ins"],
        max_ins = config["bowtie"]["max_ins"]
    output:
        bam ="alignment/{sample}.bam",
        unaligned_1 = temp("alignment/unaligned-{sample}_1.fastq"),
        unaligned_1_gz = "alignment/unaligned-{sample}_1.fastq.gz",
        unaligned_2 = temp("alignment/unaligned-{sample}_2.fastq"),
        unaligned_2_gz = "alignment/unaligned-{sample}_2.fastq.gz"
    threads: config["threads"]
    log:
       "logs/bowtie/bowtie-align-{sample}.log"
    shell: """
        (bowtie -v {params.max_mismatch} -I {params.min_ins} -X {params.max_ins} --fr --nomaqround --best -S -p {threads} --un alignment/unaligned-{wildcards.sample}.fastq {params.idx_path}/{params.basename} -1 {input.r1} -2 {input.r2} | samtools view -buh -f 0x2 - | samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {log}
        (pigz -fk alignment/unaligned-{wildcards.sample}_1.fastq alignment/unaligned-{wildcards.sample}_2.fastq) &>> {log}
        """

#the index is required use region arguments in samtools view to separate the species
rule samtools_index:
    input:
        "alignment/{sample}.bam"
    output:
        "alignment/{sample}.bam.bai"
    log:
        "logs/samtools_index/samtools_index-{sample}.log"
    shell: """
        (samtools index -b {input}) &> {log}
        """

rule bam_separate_species:
    input:
        bam = "alignment/{sample}.bam",
        bai = "alignment/{sample}.bam.bai",
        chrsizes = config["combinedgenome"]["chrsizes"]
    output:
        "alignment/{sample}_{species}only.bam"
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -bh {input.bam} $(grep {wildcards.species} {input.chrsizes} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') > {output}) &> {log}
        """

#bam must be sorted by name for bedpe. We don't do this in the bowtie step since samtools index required position-sorted bam.
rule get_fragments:
    input:
        bam = "alignment/{sample}_{species}only.bam"
    output:
        "alignment/fragments/{sample}-{species}fragments.bedpe"
    threads: config["threads"]
    log : "logs/get_fragments/get_fragments-{sample}-{species}.log"
    shell: """
        (samtools sort -n -T .{wildcards.sample}_{wildcards.species} -@ {threads} {input.bam} | bedtools bamtobed -bedpe -i stdin > {output}) &> {log}
        """

rule midpoint_coverage:
    input:
        bedpe = lambda wildcards: "alignment/fragments/" + wildcards.sample + "-" + config["combinedgenome"]["experimental_prefix"] + "fragments.bedpe" if wildcards.counttype=="counts" else "alignment/fragments/" + wildcards.sample + "-" + config["combinedgenome"]["spikein_prefix"] + "fragments.bedpe",
        chrsizes = lambda wildcards: config["genome"]["chrsizes"] if wildcards.counttype=="counts" else config["genome"]["sichrsizes"]
    params:
        prefix = lambda wildcards: config["combinedgenome"]["experimental_prefix"] if wildcards.counttype=="counts" else config["combinedgenome"]["spikein_prefix"]
    output:
        "coverage/{counttype,counts|sicounts}/{sample}-mnase-midpoint-{counttype}.bedgraph"
    log: "logs/midpoint_coverage/midpoint_coverage-{sample}-{counttype}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} width=$6-$2 {{if(width % 2 != 0){{width -= 1}}; mid=$2+width/2; print $1, mid, mid+1, $7}}' {input.bedpe} | sed -e 's/{params.prefix}//g' | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -g {input.chrsizes} -bga | sort -k1,1 -k2,2n > {output}) &> {log}
        """

# rule whole_fragment_coverage:
#     input:
#         bam = lambda wildcards: "alignment/" + wildcards.sample + "-" + config["combinedgenome"]["experimental_prefix"] + "only.bam" if wildcards.counttype=="counts" else "alignment/" + wildcards.sample + "-" + config["combinedgenome"]["spikein_prefix"] + "only.bam",
#     output:
#         "coverage/{counttype}/{sample}-mnase-wholefrag-{counttype}.bedgraph"
#     log : "logs/total_coverage/total_coverage-{sample}-{counttype}.log"
#     shell: """
#         (bedtools genomecov -ibam {input.bam} -bga -pc | sort -k1,1 -k2,2n > {output}) &> {log}
#         """

rule normalize:
    input:
        coverage = "coverage/counts/{sample}-mnase-{readtype}-counts.bedgraph",
        fragcounts = lambda wildcards: "coverage/counts/" + wildcards.sample + "-mnase-midpoint-counts.bedgraph" if wildcards.norm=="libsizenorm" else "coverage/sicounts/" + wildcards.sample + "-mnase-midpoint-sicounts.bedgraph"
    params:
        scalefactor = lambda wildcards: config["spikein-pct"] if wildcards.norm=="spikenorm" else 1
    output:
        "coverage/{norm}/{sample}-mnase-{readtype}-{norm}.bedgraph"
    wildcard_constraints:
        norm="libsizenorm|spikenorm"
    log: "logs/normalize/normalize-{sample}-{norm}-{readtype}.log"
    shell: """
        (bash scripts/libsizenorm.sh {input.fragcounts} {input.coverage} {params.scalefactor} > {output}) &> {log}
        """

rule bg_to_bw:
    input:
        bedgraph = "coverage/{norm}/{sample}-mnase-{readtype}-{norm}.bedgraph",
        chrsizes = lambda wildcards: config["genome"]["sichrsizes"] if wildcards.norm=="sicounts" else config["genome"]["chrsizes"]
    output:
        "coverage/{norm}/{sample}-mnase-{readtype}-{norm}.bw"
    log : "logs/bg_to_bw/bg_to_bw-{sample}-{readtype}-{norm}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
        """

# rule smoothed_midpoint_coverage:
#     input:
#         "coverage/bw/{sample}-mnase-midpoint-CPM.bw"
#     output:
#         "coverage/bw/{sample}-mnase-midpoint-CPM-smoothed.bw"
#     params:
#         bandwidth = config["smooth_bandwidth"]
#     log: "logs/smoothed_midpoint_coverage/smooth_midpoint_coverage-{sample}.log"
#     shell: """
#         (python scripts/smooth_midpoint_coverage.py -b {params.bandwidth} -i {input} -o {output}) &> {log}
#         """

rule deeptools_matrix:
    input:
        annotation = lambda wildcards: config["annotations"][wildcards.annotation]["path"],
        bw = "coverage/{norm}/{sample}-mnase-{readtype}-{norm}.bw"
    output:
        dtfile = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{readtype}-{norm}.mat.gz"),
        matrix = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{readtype}-{norm}.tsv"),
        matrix_gz  = "datavis/{annotation}/{norm}/{annotation}-{sample}-{readtype}-{norm}.tsv.gz"
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"] + config["annotations"][wildcards.annotation]["binsize"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"] + config["annotations"][wildcards.annotation]["binsize"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{annotation}-{sample}-{readtype}-{norm}.log"
    run:
        if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}; pigz -fk {output.matrix}) &> {log}")
        else:
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}; pigz -fk {output.matrix}) &> {log}")

rule melt_matrix:
    input:
        matrix = "datavis/{annotation}/{norm}/{annotation}-{sample}-{readtype}-{norm}.tsv.gz"
    output:
        temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{readtype}-{norm}-melted.tsv.gz")
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
    script:
        "scripts/melt_matrix.R"

rule cat_matrices:
    input:
        expand("datavis/{{annotation}}/{{norm}}/{{annotation}}-{sample}-{{readtype}}-{{norm}}-melted.tsv.gz", sample=SAMPLES)
    output:
        "datavis/{annotation}/{norm}/allsamples-{annotation}-{readtype}-{norm}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{annotation}-{readtype}-{norm}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule r_heatmaps:
    input:
        matrix = "datavis/{annotation}/{norm}/allsamples-{annotation}-{readtype}-{norm}.tsv.gz"
    output:
        heatmap_sample = "datavis/{annotation}/{norm}/mnase-{annotation}-{norm}-{status}_{condition}-v-{control}-{readtype}-heatmap-bysample.svg",
        heatmap_group = "datavis/{annotation}/{norm}/mnase-{annotation}-{norm}-{status}_{condition}-v-{control}-{readtype}-heatmap-bygroup.svg",
    params:
        samplelist = plotcorrsamples,
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
        cluster = lambda wildcards : config["annotations"][wildcards.annotation]["cluster"],
        nclust = lambda wildcards: config["annotations"][wildcards.annotation]["nclusters"],
        heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_colormap"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plot_mnase_heatmaps.R"

rule r_metagenes:
    input:
        qfrags =  "datavis/{annotation}/{norm}/allsamples-{annotation}-{readtype}-{norm}.tsv.gz"
    output:
        pmeta_group = "datavis/{annotation}/{norm}/mnase-{annotation}-{norm}-{status}_{condition}-v-{control}-{readtype}-metagene-bygroup.svg",
        pmeta_sample = "datavis/{annotation}/{norm}/mnase-{annotation}-{norm}-{status}_{condition}-v-{control}-{readtype}-metagene-bysample.svg",
        pmeta_goverlay = "datavis/{annotation}/{norm}/mnase-{annotation}-{norm}-{status}_{condition}-v-{control}-{readtype}-metagene-groupoverlay.svg",
        pmeta_soverlay = "datavis/{annotation}/{norm}/mnase-{annotation}-{norm}-{status}_{condition}-v-{control}-{readtype}-metagene-sampleoverlayall.svg",
        pmeta_soverlay_bygroup = "datavis/{annotation}/{norm}/mnase-{annotation}-{norm}-{status}_{condition}-v-{control}-{readtype}-metagene-sampleoverlaybygroup.svg"
    params:
        samplelist = plotcorrsamples,
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        trim_pct = lambda wildcards : config["annotations"][wildcards.annotation]["trim_pct"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plot_mnase_metagenes.R"

# rule build_nucwave_input:
#     input:
#        "alignment/{sample}.bam"
#     output:
#        "alignment/{sample}.bowtie"
#     log : "logs/build_nucwave_input/nucwave_input-{sample}.log"
#     shell: """
#         (samtools view {input} | awk 'BEGIN{{FS=OFS="\t"}} ($2==163) || ($2==99) {{print "+", $3, $4-1, $10}} ($2==83) || ($2==147) {{print "-", $3, $4-1, $10}}' > {output}) &> {log}
#         """

# rule nucwave:
#     input:
#         fasta = config["genome"]["fasta"],
#         alignment = "alignment/{sample}.bowtie"
#     output:
#         "nucwave/{sample}/{sample}_cut_p.wig",
#         "nucwave/{sample}/{sample}_cut_m.wig",
#         "nucwave/{sample}/{sample}_depth_complete_PE.wig",
#         "nucwave/{sample}/{sample}_PEcenter.wig",
#         "nucwave/{sample}/{sample}_depth_trimmed_PE.wig",
#         "nucwave/{sample}/{sample}_depth_wl_trimmed_PE.wig",
#         "nucwave/{sample}/{sample}.historeadsize"
#     log:
#         "logs/nucwave/nucwave-{sample}.log"
#     shell: """
#        (python scripts/nucwave_pe.py -w -o nucwave/{wildcards.sample} -g {input.fasta} -a {input.alignment} -p {wildcards.sample}) &> {log}
#        """

# rule get_fragsizes:
#     input:
#         "alignment/fragments/{sample}-fragments.bedpe"
#     output:
#         temp("alignment/fragments/.{sample}-fragsizes.tsv")
#     params:
#         group = lambda wildcards: config["samples"][wildcards.sample]["group"]
#     log : "logs/get_fragsizes/get_fragsizes-{sample}.log"
#     shell: """
#         (awk -v sample={wildcards.sample} -v group={params.group} 'BEGIN{{FS=OFS="\t"; srand()}} !/^$/ {{if (rand()<= .08) print sample, group, $6-$2}}' {input} > {output}) &> {log}
#         """

# rule cat_fragsizes:
#     input:
#         expand("alignment/fragments/.{sample}-fragsizes.tsv", sample=SAMPLES)
#     output:
#         "alignment/fragments/fragsizes.tsv.gz"
#     shell: """
#         cat {input} | pigz -f > {output}
#         """

# rule plot_fragsizes:
#     input:
#         table = "alignment/fragments/fragsizes.tsv.gz"
#     output:
#         plot = "qual_ctrl/frag-size-dist.svg"
#     script:
#         "scripts/plotfragsizedist.R"

# rule make_bigwig_for_deeptools:
#     input:
#         bg = "coverage/{sample}-mnase-midpoint-CPM.bedgraph",
#         chrsizes = config["genome"]["chrsizes"]
#     output:
#         "coverage/bw/{sample}-mnase-midpoint-CPM.bw"
#     log: "logs/make_bigwig_for_deeptools/make_bw-{sample}.log"
#     shell: """
#         (bedGraphToBigWig {input.bg} {input.chrsizes} {output}) &> {log}
#         """

# rule meta_oneoff:
#     input:
#         matrix = "datavis/allcodingTSS/allsamples-allcodingTSS.tsv.gz"
#     output:
#         metagene_overlay = "datavis/mnase-{annotation}-metagene-groupoverlay-oneoff.pdf"
#     params:
#         binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
#         upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
#         dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
#         pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
#         heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_colormap"],
#         metagene_palette = lambda wildcards : config["annotations"][wildcards.annotation]["metagene_palette"],
#         avg_heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["avg_heatmap_cmap"],
#         refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
#         ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
#     script:
#         "scripts/plotMetaOneOff.R"

# rule group_bam_for_danpos:
#     input:
#         "alignment/{sample}.bam"
#     output:
#         "danpos/{group}/{sample}.bam"
#     shell: """
#         cp {input} {output}
#     """

# rule danpos:
#     input:
#         ["danpos/" + SAMPLES[x]['group'] + "/" + x + ".bam" for x in SAMPLES]
#     output:
#         "danpos/{condition}-v-{control}/danpos_{condition}-danpos_{control}.positions.integrative.xls",
#         "danpos/{condition}-v-{control}/reference_positions.xls",
#         "danpos/{condition}-v-{control}/diff/danpos_{condition}-danpos_{control}.pois_diff.wig",
#         "danpos/{condition}-v-{control}/pooled/danpos_{condition}.Fnor.smooth.positions.ref_adjust.xls",
#         "danpos/{condition}-v-{control}/pooled/danpos_{condition}.Fnor.smooth.positions.xls",
#         "danpos/{condition}-v-{control}/pooled/danpos_{condition}.Fnor.smooth.wig",
#         "danpos/{condition}-v-{control}/pooled/danpos_{control}.Fnor.smooth.positions.ref_adjust.xls",
#         "danpos/{condition}-v-{control}/pooled/danpos_{control}.Fnor.smooth.positions.xls",
#         "danpos/{condition}-v-{control}/pooled/danpos_{control}.Fnor.smooth.wig"
#     conda:
#         "envs/danpos.yaml"
#     log: "logs/danpos/danpos-{condition}-v-{control}.log"
#     shell: """
#         module load seq/samtools/0.1.19
#         (python /groups/winston/jc459/spt5/mnase-seq/scripts/danpos-2.2.2/danpos.py dpos danpos/{wildcards.condition}/:danpos/{wildcards.control}/ -m 1 -o danpos/{wildcards.condition}-v-{wildcards.control}) &> {log}
#         module unload seq/samtools/0.1.19
#         """

# ##for now I use default parameters except for paired end
