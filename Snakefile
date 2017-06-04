#!/usr/bin/env python

configfile: "config.yaml"

SAMPLES = config["barcodes"]

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
            melt_matrix,
            gzip_deeptools_table
#requirements
# seq/ea-utils/
# UNLOAD seq/cutadapt/1.11
# use pip-installed cutadapt
# seq/bowtie/1.1.1
# seq/samtools/1.3
# pip install numpy PyWavelets
# pip install deeptools
# UNLOAD deeptools, pysam
# UNLOAD stats/R/3.2.1, load stats/R/3.3.1
rule all:
    input:
        "qual_ctrl/raw",
        #"qual_ctrl/coverage.tsv",
        #"qual_ctrl/coverage.png",
        "qual_ctrl/frag-size-dist.png",
        "qual_ctrl/seq-depth-dist.png",
        expand("nucwave/{sample}/{sample}_depth_wl_trimmed_PE.wig", sample=SAMPLES),
        expand("coverage/{sample}-midpoint-CPM.bedgraph", sample=SAMPLES),
        "correlations/pca-scree.png",
        expand("heatmaps/{annotation}/{annotation}-{sample}-heatmap.png", annotation=config["annotations"], sample=SAMPLES),
        expand("alignment/unaligned-{sample}_2.fastq.gz", sample=SAMPLES),
        expand("metagene/{annotation}/{annotation}-{sample}-metagene.png", annotation = config["annotations"], sample=SAMPLES),
        #expand("metagene/{annotation}/{annotation}-allsamples.png", annotation = config["annotations"])

rule make_barcode_file:
    output:
        "barcodes.tsv"
    params:
        bc = config["barcodes"]
    run:
        with open(output[0], "w") as out:
            for x in params.bc:
                out.write(('\t'.join((x, params.bc[x]))+'\n'))

#fastQC on raw sequencing data
rule fastqc_raw:
    input:
       r1 = config["fastq"]["r1"],
       r2 = config["fastq"]["r2"]
    output:
       "qual_ctrl/raw"
    threads : config["threads"]
    log : "logs/fastqc/fastqc-raw.log"
    shell: """
        mkdir -p {output}
        (fastqc -o {output} --noextract -t {threads} {input.r1} {input.r2}) &> {log}
        """

rule demultiplex:
    input:
        r1 = config["fastq"]["r1"],
        r2 = config["fastq"]["r2"],
        barcodes = "barcodes.tsv"
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
#    note: the minimum length requirement (trimmed read >= 5nt) is important to sanitize the output for bowtie 1
rule cutadapt:
    input:
        r1 = "fastq/{sample}.r1.fastq.gz",
        r2 = "fastq/{sample}.r2.fastq.gz"
    output:
        r1 = temp("fastq/trimmed/{sample}-trim.r1.fastq"),
        r2 = temp("fastq/trimmed/{sample}-trim.r2.fastq")
    params:
        qual_cutoff = config["cutadapt"]["qual_cutoff"],
        adapter = lambda wildcards : SAMPLES[wildcards.sample]+"T"
    log:
        "logs/cutadapt/cutadapt-{sample}.log"
    shell: """
        (cutadapt -u 1 -G ^{params.adapter} -q {params.qual_cutoff} --minimum-length 5 -o {output.r1} -p {output.r2} {input.r1} {input.r2}) &> {log}
        """

#fastQC on demultiplexed and trimmed reads
rule fastqc_processed:
    input:
        r1 = "fastq/trimmed/{sample}-trim.r1.fastq",
        r2 = "fastq/trimmed/{sample}-trim.r2.fastq"
    output:
        "qual_ctrl/trim-{sample}/{sample}-trim.r1_fastqc.html",
        "qual_ctrl/trim-{sample}/{sample}-trim.r2_fastqc.html"
    threads : config["threads"]
    log : "logs/fastqc/fastqc-trim-{sample}.log"
    shell: """
        mkdir -p {output}
        (fastqc -o qual_ctrl/trim-{wildcards.sample} --noextract -t {threads} {input.r1} {input.r2}) &> {log}
        """

#build bowtie index for genome
basename = config["genome"]["name"]

rule bowtie_build:
    input:
        fasta = config["genome"]["fasta"]
    output:
        expand("../genome/bowtie1_indexes/" + basename + ".{n}.ebwt", n=[1,2,3,4]),
        expand("../genome/bowtie1_indexes/" + basename + ".rev.{n}.ebwt", n=[1,2])
    params:
        outbase = "../genome/bowtie1_indexes/" + basename
    log:
        "logs/bowtie/bowtie-build.log"    
    shell: """
        (bowtie-build {input.fasta} {params.outbase}) &> {log}
        """

#align with Bowtie 1
#in Christine's paper, Burak uses -m 10 --best
rule bowtie:
    input:
        expand("../genome/bowtie1_indexes/" + basename + ".{n}.ebwt", n=[1,2,3,4]),
        expand("../genome/bowtie1_indexes/" + basename + ".rev.{n}.ebwt", n=[1,2]),
        r1 = "fastq/trimmed/{sample}-trim.r1.fastq",
        r2 = "fastq/trimmed/{sample}-trim.r2.fastq"
    output:
        bam ="alignment/{sample}.bam",
        unaligned_1 = temp("alignment/unaligned-{sample}_1.fastq"),
        unaligned_2 = temp("alignment/unaligned-{sample}_2.fastq")
    threads: config["threads"]
    params:
        outbase = "../genome/bowtie1_indexes/" + basename,
        max_mismatch = config["bowtie"]["max_mismatch"],
        min_ins = config["bowtie"]["min_ins"],
        max_ins = config["bowtie"]["max_ins"]
    log:
       "logs/bowtie/bowtie-align-{sample}.log"
    shell: """
        (bowtie -v {params.max_mismatch} -I {params.min_ins} -X {params.max_ins} --fr --nomaqround --best -S -p {threads} --un alignment/unaligned-{wildcards.sample}.fastq {params.outbase} -1 {input.r1} -2 {input.r2} | samtools view -buh -f 0x2 - | samtools sort -T {wildcards.sample} -@ {threads} -o {output.bam} -) &> {log}
        """

rule gzip_loose_fastq:
    input:
        "alignment/{sample}.bam",
        "qual_ctrl/trim-{sample}/{sample}-trim.r1_fastqc.html",
        "qual_ctrl/trim-{sample}/{sample}-trim.r2_fastqc.html",
        trim_r1 = "fastq/trimmed/{sample}-trim.r1.fastq",
        trim_r2 = "fastq/trimmed/{sample}-trim.r2.fastq",
        unaligned_r1 = "alignment/unaligned-{sample}_1.fastq",
        unaligned_r2 = "alignment/unaligned-{sample}_2.fastq"
    output:
        trim_r1 = "fastq/trimmed/{sample}-trim.r1.fastq.gz",
        trim_r2 = "fastq/trimmed/{sample}-trim.r2.fastq.gz",
        unaligned_r1 = "alignment/unaligned-{sample}_1.fastq.gz",
        unaligned_r2 = "alignment/unaligned-{sample}_2.fastq.gz"
    threads : config["threads"]
    shell: """
        pigz -fk {input.trim_r1}
        pigz -fk {input.trim_r2}
        pigz -fk {input.unaligned_r1}
        pigz -fk {input.unaligned_r2}
        """
        

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

rule deeptools_plotcoverage:
    input:
        bam = expand("alignment/{sample}.bam", sample=SAMPLES),
        index = expand("alignment/{sample}.bam.bai", sample=SAMPLES)
    output:
        table= "qual_ctrl/coverage.tsv",
        fig="qual_ctrl/coverage.png"
    params:
        min_ins = config["bowtie"]["min_ins"],
        max_ins = config["bowtie"]["max_ins"]
    log: "logs/deeptools/deeptools_plotcoverage.log"
    threads: config["threads"]
    shell: """
        (plotCoverage --outRawCounts {output.table} -p {threads} -v --extendReads --minFragmentLength {params.min_ins} --maxFragmentLength {params.max_ins} -b {input.bam} -o {output.fig}) &> {log}
        """   

rule build_nucwave_input:
    input:
       "alignment/{sample}.bam" 
    output:
       "alignment/{sample}.bowtie"
    log : "logs/build_nucwave_input/nucwave_input-{sample}.log"
    shell: """
        (samtools view {input} | awk 'BEGIN{{FS=OFS="\t"}} ($2==163) || ($2==99) {{print "+", $3, $4-1, $10}} ($2==83) || ($2==147) {{print "-", $3, $4-1, $10}}' > {output}) &> {log}
        """

rule nucwave:
    input:
        fasta = config["genome"]["fasta"],
        alignment = "alignment/{sample}.bowtie"
    output:
        "nucwave/{sample}/{sample}_cut_p.wig",
        "nucwave/{sample}/{sample}_cut_m.wig",
        "nucwave/{sample}/{sample}_depth_complete_PE.wig",
        "nucwave/{sample}/{sample}_PEcenter.wig",
        "nucwave/{sample}/{sample}_depth_trimmed_PE.wig",
        "nucwave/{sample}/{sample}_depth_wl_trimmed_PE.wig",
        "nucwave/{sample}/{sample}.historeadsize"
    log:
        "logs/nucwave/nucwave-{sample}.log"
    shell: """
       (python scripts/nucwave_pe.py -w -o nucwave/{wildcards.sample} -g {input.fasta} -a {input.alignment} -p {wildcards.sample}) &> {log}
       """

rule get_fragments:
    input:
        bam = "alignment/{sample}.bam"
    output:
        "alignment/fragments/{sample}-fragments.bedpe"
    threads: config["threads"]
    log : "logs/get_fragments/get_fragments-{sample}.log"
    shell: """
        (samtools sort -n -@ {threads} {input.bam} | bedtools bamtobed -bedpe -i stdin > {output}) &> {log}
        """

rule get_fragsizes:
    input:
        "alignment/fragments/{sample}-fragments.bedpe"
    output:
        temp("alignment/fragments/.{sample}-fragsizes.tsv")
    log : "logs/get_fragsizes/get_fragsizes-{sample}.log"
    shell: """
        (awk -v sample={wildcards.sample} 'BEGIN{{FS=OFS="\t"; srand()}} !/^$/ {{if (rand()<= .08) print sample , $6-$2}}' {input} > {output}) &> {log}
        """

rule cat_fragsizes:
    input:
        expand("alignment/fragments/.{sample}-fragsizes.tsv", sample=SAMPLES)
    output:
        "alignment/fragments/fragsizes.tsv"
    shell: """
        cat {input} > {output}
        """

rule plot_fragsizes:
    input:
        table = "alignment/fragments/fragsizes.tsv"
    output:
        plot = "qual_ctrl/frag-size-dist.png"
    script:
        "scripts/plotfragsizedist.R"

rule make_fragments_table:
    input:
        expand("alignment/fragments/{sample}-fragments.bedpe", sample=SAMPLES)
    output:
        "qual_ctrl/fragment_counts.txt"
    shell: """
        wc -l {input} > {output}
        """
#note: to retrieve number of fragments in a sample, | grep {sample} qual_ctrl/fragment_counts.txt | cut -f3 -d ' '

rule get_perbase_depth:
    input:
        bam = "alignment/{sample}.bam",
        chrsizes = config["genome"]["chrsizes"]
    output:
        temp("qual_ctrl/.{sample}-depth.tsv")
    log: "logs/get_perbase_depth/perbase_depth-{sample}.log"
    shell: """
        (bedtools genomecov -ibam {input.bam} -g {input.chrsizes} -d -pc | awk -v sample={wildcards.sample} 'BEGIN{{FS=OFS="\t"; srand()}} !/^$/ {{if (rand()<= .08) print sample, $3}}' > {output}) &> {log}
        """

rule cat_perbase_depth:
    input:
        expand("qual_ctrl/.{sample}-depth.tsv", sample=SAMPLES)
    output:
        "qual_ctrl/perbasedepth.tsv.gz"
    #params:
    #    labels = list(SAMPLES.keys())
    shell: """
        cat {input} | gzip -f > {output}  
        """

rule plot_depth:
    input:
        "qual_ctrl/perbasedepth.tsv.gz"
    output:
        "qual_ctrl/seq-depth-dist.png"
    script:
        "scripts/plotCoverage.R"


rule midpoint_coverage:
    input:
        bedpe = "alignment/fragments/{sample}-fragments.bedpe",
        chrsizes = config["genome"]["chrsizes"]
    output:
        "coverage/{sample}-midpoint-counts.bedgraph"
    log: "logs/midpoint_coverage/midpoint_coverage-{sample}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} width=$6-$2 {{if(width % 2 != 0){{width -= 1}}; mid=$2+width/2; print $1, mid, mid+1, $7}}' {input.bedpe} | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -g {input.chrsizes} -bga | sortBed -i stdin > {output}) &> {log}
        """

rule total_coverage:
    input:
        bam = "alignment/{sample}.bam",
        chrsizes = config["genome"]["chrsizes"]
    output:
        "coverage/{sample}-total-counts.bedgraph"
    log : "logs/total_coverage/total_coverage-{sample}.log"
    shell: """
        (bedtools genomecov -ibam {input.bam} -g {input.chrsizes} -bga -pc | sortBed -i stdin > {output}) &> {log}
        """

rule seq_depth_norm:
    input:
        midpoint = "coverage/{sample}-midpoint-counts.bedgraph",
        total = "coverage/{sample}-total-counts.bedgraph",
        nfrags = "qual_ctrl/fragment_counts.txt"
    output:
        midpoint = "coverage/{sample}-midpoint-CPM.bedgraph",
        total = "coverage/{sample}-total-CPM.bedgraph"
    log : "logs/seq_depth_norm/depthnorm-{sample}.log"
    shell: """
        norm=$(grep /{wildcards.sample} {input.nfrags} | awk '{{print $1}}')
        (awk -v anorm=$norm 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4*1000000/anorm}}' {input.midpoint} > {output.midpoint}) &>> {log}
        (awk -v anorm=$norm 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4*1000000/anorm}}' {input.total} > {output.total}) &>> {log}
        """

rule make_window_file:
    input:
        chrsizes = config["genome"]["chrsizes"]
    output:
        temp("genome/windows.bed")
    params:
        wsize = config["corr-binsize"]
    log: "logs/bedtools_make_window_file.log"
    shell: """
        (bedtools makewindows -g {input.chrsizes} -w {params.wsize} | sortBed -i stdin > {output}) &> {log}
        """

rule map_to_windows:
    input:
        bed = "genome/windows.bed",
        bedgraph = "coverage/{sample}-midpoint-CPM.bedgraph"
    output:
        temp("coverage/.{sample}-maptowindow.tsv")
    log : "logs/bedtools_map_to_windows/bedtools_map-{sample}.log"
    shell: """
        (bedtools map -c 4 -o sum -a {input.bed} -b {input.bedgraph} | cut -f4 > {output}) &> {log}
        """

rule cat_windows:
    input:
        values = expand("coverage/.{sample}-maptowindow.tsv", sample=SAMPLES),
        coord = "genome/windows.bed"
    output:
        "correlations/midpoint-CPM-windows.tsv"
    params:
        labels = list(config["barcodes"].keys()),
    shell: """
        echo -e "chr\tstart\tend\t{params.labels}\n$(paste {input.coord} {input.values})" > {output}
        """

rule plot_correlations:
    input:
        "correlations/midpoint-CPM-windows.tsv"
    output:
        scatter = "correlations/pairwise-scatterplots.png",
        dists_cluster = "correlations/sample-dists-clustered.png",
        dists_nocluster = "correlations/sample-dists-unclustered.png",
        pca = "correlations/pca.png",
        scree = "correlations/pca-scree.png"
    script:
        "scripts/plotcorrelations.R"
                
rule make_bigwig_for_deeptools:
    input:
        bg = "coverage/{sample}-midpoint-CPM.bedgraph",
        chrsizes = config["genome"]["chrsizes"]
    output:
        "coverage/bw/{sample}-midpoint-CPM.bw"
    log: "logs/make_bigwig_for_deeptools/make_bw-{sample}.log"
    shell: """
        (bedGraphToBigWig {input.bg} {input.chrsizes} {output}) &> {log}
        """

rule deeptools_matrix:
    input:
        annotation = lambda wildcards: config["annotations"][wildcards.annotation]["path"],
        bw = "coverage/bw/{sample}-midpoint-CPM.bw"
    output:
        dtfile = "heatmaps/{annotation}/{annotation}-{sample}.mat.gz",
        matrix = temp("heatmaps/{annotation}/{annotation}-{sample}.tsv")
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{annotation}-{sample}.log"
    shell: """
        (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}
        """
    #shell: """
    #    (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --sortRegions {params.sort} --averageTypeBins {params.binstat} -p {threads}) &> {log}
    #    """

rule gzip_deeptools_table:
    input:
        "heatmaps/{annotation}/{annotation}-{sample}.tsv"
    output:
        "heatmaps/{annotation}/{annotation}-{sample}.tsv.gz"
    shell: """
        pigz -f {input}
        """

rule r_plotHeatmap:
    input:
        matrix = "heatmaps/{annotation}/{annotation}-{sample}.tsv.gz"
    output:
        "heatmaps/{annotation}/{annotation}-{sample}-heatmap.png"
    params:
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        figunits = lambda wildcards : config["annotations"][wildcards.annotation]["figunits"],
        figwidth = lambda wildcards : config["annotations"][wildcards.annotation]["figwidth"],
        figheight = lambda wildcards : config["annotations"][wildcards.annotation]["figheight"],
        cmap = lambda wildcards : config["annotations"][wildcards.annotation]["colormap"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plotHeatmap.R"
 
rule melt_matrix:
    input:
        matrix = "heatmaps/{annotation}/{annotation}-{sample}.tsv.gz"
    output:
        "metagene/{annotation}/{annotation}-{sample}-melted.tsv.gz"
    params:
        name = lambda wildcards : wildcards.sample,
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
    script:
        "scripts/melt_matrix.R"
    
rule cat_matrices:
    input:
        expand("metagene/{{annotation}}/{{annotation}}-{sample}-melted.tsv.gz", sample=SAMPLES)
    output:
        "metagene/{annotation}/{annotation}-allsamples.tsv.gz"
    shell: """
        cat {input} > {output}
        """

rule plot_indiv_metagene:
    input:
        "metagene/{annotation}/{annotation}-{sample}-melted.tsv.gz"
    output:
        "metagene/{annotation}/{annotation}-{sample}-metagene.png"
    params:
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"]
    script:
        "scripts/plotMetagene2.R"

rule plot_combined_metagene:
    input:
        "metagene/{annotation}/{annotation}-allsamples.tsv.gz"
    output:
        "metagene/{annotation}/{annotation}-allsamples.png"
    params:
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"]
    script:
        "scripts/plotMetagene.R"

