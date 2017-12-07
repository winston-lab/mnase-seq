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
            build_nucwave_input,
            get_fragsizes,
            cat_fragsizes,
            make_fragments_table,
            make_window_file,
            cat_windows,
            cat_perbase_depth,
            seq_depth_norm,
            make_bigwig_for_deeptools,
            group_bam_for_danpos,
            dpos_wig_to_bigwig,
            dpos_gzip_deeptools_table,
            dpos_melt_matrix,
            dpos_cat_matrices

rule all:
    input:
        #fastqc
        "qual_ctrl/fastqc/raw",
        # expand("qual_ctrl/fastqc/{sample}/{sample}-clean.{read}_fastqc.html", sample=SAMPLES, read=["r1","r2"]),
        #demultiplex
        expand("fastq/{sample}.{read}.fastq.gz", sample=SAMPLES, read=["r1","r2"]),
        #alignment
        # expand("alignment/{sample}.bam", sample=SAMPLES),
        # expand("alignment/unaligned-{sample}_{r}.fastq.gz", sample=SAMPLES, r=[1,2])

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
#TODO: check this works for barcodes of different lengths
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
rule cutadapt:
    input:
        r1 = "fastq/{sample}.r1.fastq.gz",
        r2 = "fastq/{sample}.r2.fastq.gz"
    output:
        r1 = "fastq/cleaned/{sample}-clean.r1.fastq.gz",
        r2 = "fastq/cleaned/{sample}-clean.r2.fastq.gz"
    params:
        qual_cutoff = config["cutadapt"]["qual_cutoff"],
        adapter = lambda wildcards : SAMPLES[wildcards.sample]["barcode"]+"T"
    log:
        "logs/cutadapt/cutadapt-{sample}.log"
    shell: """
        (cutadapt -u 1 -G ^{params.adapter} -q {params.qual_cutoff} --minimum-length 5 -o {output.r1} -p {output.r2} {input.r1} {input.r2}) &> {log}
        """

#fastQC on demultiplexed and cleaned reads
rule fastqc_processed:
    input:
        r1 = "fastq/cleaned/{sample}-clean.r1.fastq.gz",
        r2 = "fastq/cleaned/{sample}-clean.r2.fastq.gz"
    output:
        "qual_ctrl/fastqc/{sample}/{sample}-clean.r1_fastqc.html",
        "qual_ctrl/fastqc/{sample}/{sample}-clean.r2_fastqc.html"
    threads : config["threads"]
    log : "logs/fastqc/fastqc-trim-{sample}.log"
    shell: """
        mkdir -p {output}
        (fastqc -o qual_ctrl/trim-{wildcards.sample} --noextract -t {threads} {input.r1} {input.r2}) &> {log}
        """

rule bowtie_build:
    input:
        fasta = config["genome"]["fasta"]
    output:
        expand(config["bowtie"]["index-path"] + "/{{basename}}.{num}.ebwt", num=[1,2,3,4]),
        expand(config["bowtie"]["index-path"] + "/{{basename}}.rev.{num}.ebwt", num=[1,2]),
    params:
        idx_path = config["bowtie"]["index-path"],
        prefix = config["combinedgenome"]["experimental_prefix"]
    log:
        "logs/bowtie/bowtie-build.log"
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
        expand(config["bowtie"]["index-path"] + "/" + config["combinedgenome"]["name"] + "rev.{num}.ebwt", num=[1,2]) if sisamples else expand(config["bowtie"]["index-path"] + "/" + config["genome"]["name"] + "rev.{num}.ebwt", num=[1,2]),
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
        unaligned_2 = temp("alignment/unaligned-{sample}_2.fastq")
    threads: config["threads"]
    log:
       "logs/bowtie/bowtie-align-{sample}.log"
    shell: """
        (bowtie -v {params.max_mismatch} -I {params.min_ins} -X {params.max_ins} --fr --nomaqround --best -S -p {threads} --un alignment/unaligned-{wildcards.sample}.fastq {params.idx_path}/{params.basename} -1 {input.r1} -2 {input.r2} | samtools view -buh -f 0x2 - | samtools sort -T {wildcards.sample} -@ {threads} -o {output.bam} -) &> {log}
        """

rule gzip_loose_fastq:
    input:
        "alignment/{sample}.bam",
        unaligned_r1 = "alignment/unaligned-{sample}_1.fastq",
        unaligned_r2 = "alignment/unaligned-{sample}_2.fastq"
    output:
        unaligned_r1 = "alignment/unaligned-{sample}_1.fastq.gz",
        unaligned_r2 = "alignment/unaligned-{sample}_2.fastq.gz"
    threads : config["threads"]
    shell: """
        pigz -f {input.unaligned_r1}
        pigz -f {input.unaligned_r2}
        """

# rule get_fragments:
#     input:
#         bam = "alignment/{sample}.bam"
#     output:
#         "alignment/fragments/{sample}-fragments.bedpe"
#     threads: config["threads"]
#     log : "logs/get_fragments/get_fragments-{sample}.log"
#     shell: """
#         (bedtools bamtobed -bedpe -i {input.bam} > {output}) &> {log}
#         """
#         # (samtools sort -n -@ {threads} {input.bam} | bedtools bamtobed -bedpe -i stdin > {output}) &> {log}

# rule midpoint_coverage:
#     input:
#         bedpe = "alignment/fragments/{sample}-fragments.bedpe",
#         chrsizes = config["genome"]["chrsizes"]
#     output:
#         "coverage/{sample}-mnase-midpoint-counts.bedgraph"
#     log: "logs/midpoint_coverage/midpoint_coverage-{sample}.log"
#     shell: """
#         (awk 'BEGIN{{FS=OFS="\t"}} width=$6-$2 {{if(width % 2 != 0){{width -= 1}}; mid=$2+width/2; print $1, mid, mid+1, $7}}' {input.bedpe} | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -g {input.chrsizes} -bga | sortBed -i stdin > {output}) &> {log}
#         """

# rule total_coverage:
#     input:
#         bam = "alignment/{sample}.bam",
#         chrsizes = config["genome"]["chrsizes"]
#     output:
#         "coverage/{sample}-mnase-total-counts.bedgraph"
#     log : "logs/total_coverage/total_coverage-{sample}.log"
#     shell: """
#         (bedtools genomecov -ibam {input.bam} -g {input.chrsizes} -bga -pc | sortBed -i stdin > {output}) &> {log}
#         """


# rule samtools_index:
#     input:
#         "alignment/{sample}.bam"
#     output:
#         "alignment/{sample}.bam.bai"
#     log:
#         "logs/samtools_index/samtools_index-{sample}.log"
#     shell: """
#         (samtools index -b {input}) &> {log}
#         """

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

# rule get_perbase_depth:
#     input:
#         bam = "alignment/{sample}.bam",
#         chrsizes = config["genome"]["chrsizes"]
#     output:
#         temp("qual_ctrl/.{sample}-depth.tsv")
#     log: "logs/get_perbase_depth/perbase_depth-{sample}.log"
#     shell: """
#         (bedtools genomecov -ibam {input.bam} -g {input.chrsizes} -d -pc | awk -v sample={wildcards.sample} 'BEGIN{{FS=OFS="\t"; srand()}} !/^$/ {{if (rand()<= .08) print sample, $3}}' > {output}) &> {log}
#         """

# rule cat_perbase_depth:
#     input:
#         expand("qual_ctrl/.{sample}-depth.tsv", sample=SAMPLES)
#     output:
#         "qual_ctrl/perbasedepth.tsv.gz"
#     shell: """
#         cat {input} | gzip -f > {output}
#         """

# rule plot_depth:
#     input:
#         "qual_ctrl/perbasedepth.tsv.gz"
#     output:
#         "qual_ctrl/seq-depth-dist.svg"
#     script:
#         "scripts/plotCoverage.R"

# rule seq_depth_norm:
#     input:
#         midpoint = "coverage/{sample}-mnase-midpoint-counts.bedgraph",
#         total = "coverage/{sample}-mnase-total-counts.bedgraph",
#         nfrags = "qual_ctrl/fragment_counts.txt"
#     output:
#         midpoint = "coverage/{sample}-mnase-midpoint-CPM.bedgraph",
#         total = "coverage/{sample}-mnase-total-CPM.bedgraph"
#     log : "logs/seq_depth_norm/depthnorm-{sample}.log"
#     shell: """
#         norm=$(grep /{wildcards.sample} {input.nfrags} | awk '{{print $1}}')
#         (awk -v anorm=$norm 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4*1000000/anorm}}' {input.midpoint} > {output.midpoint}) &>> {log}
#         (awk -v anorm=$norm 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4*1000000/anorm}}' {input.total} > {output.total}) &>> {log}
#         """

# rule make_window_file:
#     input:
#         chrsizes = config["genome"]["chrsizes"]
#     output:
#         temp("genome/windows.bed")
#     params:
#         wsize = config["corr-binsize"]
#     log: "logs/bedtools_make_window_file.log"
#     shell: """
#         (bedtools makewindows -g {input.chrsizes} -w {params.wsize} | sortBed -i stdin > {output}) &> {log}
#         """

# rule map_to_windows:
#     input:
#         bed = "genome/windows.bed",
#         bedgraph = "coverage/{sample}-midpoint-CPM.bedgraph"
#     output:
#         temp("coverage/.{sample}-maptowindow.tsv")
#     log : "logs/bedtools_map_to_windows/bedtools_map-{sample}.log"
#     shell: """
#         (bedtools map -c 4 -o sum -a {input.bed} -b {input.bedgraph} | cut -f4 > {output}) &> {log}
#         """

# rule cat_windows:
#     input:
#         values = expand("coverage/.{sample}-maptowindow.tsv", sample=SAMPLES),
#         coord = "genome/windows.bed"
#     output:
#         "correlations/midpoint-CPM-windows.tsv"
#     params:
#         labels = list(SAMPLES.keys()),
#     shell: """
#         echo -e "chr\tstart\tend\t{params.labels}\n$(paste {input.coord} {input.values})" > {output}
#         """

# rule plot_correlations:
#     input:
#         "correlations/midpoint-CPM-windows.tsv"
#     output:
#         scatter = "correlations/pairwise-scatterplots.svg",
#         dists_cluster = "correlations/sample-dists-clustered.svg",
#         dists_nocluster = "correlations/sample-dists-unclustered.svg",
#         pca = "correlations/pca.svg",
#         scree = "correlations/pca-scree.svg"
#     script:
#         "scripts/plotcorrelations.R"

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

# rule deeptools_matrix:
#     input:
#         annotation = lambda wildcards: config["annotations"][wildcards.annotation]["path"],
#         #bw = "coverage/bw/{sample}-mnase-midpoint-CPM.bw"
#         bw = "coverage/bw/{sample}-mnase-midpoint-CPM-smoothed.bw"
#     output:
#         dtfile = temp("datavis/{annotation}/{annotation}-{sample}.mat.gz"),
#         matrix = temp("datavis/{annotation}/{annotation}-{sample}.tsv")
#     params:
#         refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
#         upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
#         dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
#         binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
#         sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
#         sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
#         binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
#     threads : config["threads"]
#     log: "logs/deeptools/computeMatrix-{annotation}-{sample}.log"
#     run:
#         if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
#             shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")
#         else:
#             shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")

# rule gzip_deeptools_matrix:
#     input:
#         matrix = "datavis/{annotation}/{annotation}-{sample}.tsv"
#     output:
#         "datavis/{annotation}/{annotation}-{sample}.tsv.gz"
#     shell: """
#         pigz -f {input}
#         """

# rule melt_matrix:
#     input:
#         matrix = "datavis/{annotation}/{annotation}-{sample}.tsv.gz"
#     output:
#         temp("datavis/{annotation}/{annotation}-{sample}-melted.tsv.gz")
#     params:
#         group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
#         binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
#         upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
#         dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"]
#     script:
#         "scripts/melt_matrix.R"

# rule cat_matrices:
#     input:
#         expand("datavis/{{annotation}}/{{annotation}}-{sample}-melted.tsv.gz", sample=SAMPLES)
#     output:
#         "datavis/{annotation}/allsamples-{annotation}.tsv.gz"
#     shell: """
#         cat {input} > {output}
#         """

# rule r_datavis:
#     input:
#         matrix = "datavis/{annotation}/allsamples-{annotation}.tsv.gz"
#     output:
#         heatmap_sample = "datavis/{annotation}/mnase-{annotation}-heatmap-bysample.svg",
#         heatmap_group = "datavis/{annotation}/mnase-{annotation}-heatmap-bygroup.svg",
#         metagene_sample = "datavis/{annotation}/mnase-{annotation}-metagene-bysample.svg",
#         metagene_sample_overlay = "datavis/{annotation}/mnase-{annotation}-metagene-sampleolaybygroup.svg",
#         metagene_sample_overlay_all = "datavis/{annotation}/mnase-{annotation}-metagene-sampleolayall.svg",
#         metagene_group = "datavis/{annotation}/mnase-{annotation}-metagene-bygroup.svg",
#         metagene_overlay = "datavis/{annotation}/mnase-{annotation}-metagene-groupolay.svg",
#         metaheatmap_sample = "datavis/{annotation}/mnase-{annotation}-metaheatmap-bysample.svg",
#         metaheatmap_group = "datavis/{annotation}/mnase-{annotation}-metaheatmap-bygroup.svg"
#     params:
#         upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
#         dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
#         pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
#         trim_pct = lambda wildcards : config["annotations"][wildcards.annotation]["trim_pct"],
#         heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_colormap"],
#         metagene_palette = lambda wildcards : config["annotations"][wildcards.annotation]["metagene_palette"],
#         avg_heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["avg_heatmap_cmap"],
#         refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
#         ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
#     script:
#         "scripts/plotHeatmapsMeta.R"

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

# rule dpos_wig_to_bigwig:
#     input:
#         wig = "danpos/{condition}-v-{control}/diff/danpos_{condition}-danpos_{control}.pois_diff.wig",
#         chrsizes = config["genome"]["chrsizes"]
#     output:
#         "danpos/{condition}-v-{control}/diff/danpos_{condition}-danpos_{control}.pois_diff.bw",
#     log: "logs/dpos_wig_to_bigwig/dpos_wig_to_bigwig-{condition}-v-{control}.log"
#     shell: """
#         (wigToBigWig {input.wig} {input.chrsizes} {output} ) &> {log}
#         """

# rule dpos_deeptools_matrix:
#     input:
#         annotation = lambda wildcards: config["annotations"][wildcards.annotation]["path"],
#         bw = "danpos/{condition}-v-{control}/diff/danpos_{condition}-danpos_{control}.pois_diff.bw",
#     output:
#         dtfile = temp("datavis/{annotation}/dpos/dpos-{annotation}-{condition}-v-{control}.mat.gz"),
#         matrix = temp("datavis/{annotation}/dpos/dpos-{annotation}-{condition}-v-{control}.tsv")
#     params:
#         refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
#         upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
#         dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
#         binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
#         sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
#         sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
#         binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
#     threads : config["threads"]
#     log: "logs/deeptools/dpos_computeMatrix-{annotation}-{condition}-v-{control}.log"
#     shell: """
#         (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}
#         """
#     #shell: """
#     #    (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}
#     #    """

# rule dpos_gzip_deeptools_table:
#     input:
#         tsv = "datavis/{annotation}/dpos/dpos-{annotation}-{condition}-v-{control}.tsv",
#         mat = "datavis/{annotation}/dpos/dpos-{annotation}-{condition}-v-{control}.mat.gz"
#     output:
#         "datavis/{annotation}/dpos/dpos-{annotation}-{condition}-v-{control}-t.tsv.gz"
#     shell: """
#         pigz -fc {input.tsv} > {output}
#         rm {input.mat}
#         """

# rule dpos_melt_matrix:
#     input:
#         matrix = "datavis/{annotation}/dpos/dpos-{annotation}-{condition}-v-{control}-t.tsv.gz"
#     output:
#         temp("datavis/{annotation}/dpos/dpos-{annotation}-{condition}-v-{control}-melted.tsv.gz")
#     params:
#         controlgroup = lambda wildcards : wildcards.control,
#         conditiongroup = lambda wildcards : wildcards.condition,
#         binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
#         upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
#         dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"]
#     script:
#         "scripts/melt_lfc_matrix.R"

# rule dpos_cat_matrices:
#     input:
#         expand("datavis/{{annotation}}/dpos/dpos-{{annotation}}-{condition}-v-{control}-melted.tsv.gz", condition = conditiongroups, control = controlgroups)
#     output:
#         "datavis/{annotation}/dpos/dpos-allconditions-{annotation}.tsv.gz"
#     shell: """
#         cat {input} > {output}
#         """
