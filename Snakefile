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
            midpoint_coverage,
            total_coverage,
            seq_depth_norm,
            make_window_file,
            map_to_windows,
            cat_windows
#requirements
# seq/ea-utils/
# UNLOAD seq/cutadapt/1.11
# use pip-installed cutadapt
# seq/bowtie/1.1.1
# seq/samtools/1.3
# pip install numpy PyWavelets
# pip install deeptools
# UNLOAD deeptools, pysam

rule all:
    input:
        "qual_ctrl/raw",
        expand("qual_ctrl/trim-{sample}", sample=SAMPLES), 
        "qual_ctrl/coverage.tsv",
        "qual_ctrl/coverage.png",
        "qual_ctrl/frag-size-dist.png",
        expand("nucwave/{sample}/{sample}_depth_wl_trimmed_PE.wig", sample=SAMPLES),
        expand("coverage/{sample}-midpoint-CPM.bedgraph", sample=SAMPLES),

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
        mkdir {output}
        (fastqc -o {output} --noextract -t {threads} {input.r1} {input.r2} ) &> {log}
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
        r1 = "fastq/trimmed/{sample}-trim.r1.fastq",
        r2 = "fastq/trimmed/{sample}-trim.r2.fastq"
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
        "qual_ctrl/trim-{sample}"
    threads : config["threads"]
    log : "logs/fastqc/fastqc-trim-{sample}.log"
    shell: """
        mkdir {output}
        (fastqc -o {output} --noextract -t {threads} {input.r1} {input.r2}) &> {log}
        """

#build bowtie index for genome
basename = config["genome"]["name"]

rule bowtie_build:
    input:
        fasta = config["genome"]["fasta"]
    output:
        expand("genome/" + basename + ".{n}.ebwt", n=[1,2,3,4]),
        expand("genome/" + basename + ".rev.{n}.ebwt", n=[1,2])
    params:
        outbase = "genome/" + basename
    log:
        "logs/bowtie/bowtie-build.log"    
    shell: """
        (bowtie-build {input.fasta} {params.outbase}) &> {log}
        """

#align with Bowtie 1
#in Christine's paper, Burak uses -m 10 --best
rule bowtie:
    input:
        expand("genome/" + basename + ".{n}.ebwt", n=[1,2,3,4]),
        expand("genome/" + basename + ".rev.{n}.ebwt", n=[1,2]),
        r1 = "fastq/trimmed/{sample}-trim.r1.fastq",
        r2 = "fastq/trimmed/{sample}-trim.r2.fastq"
    output:
        bam ="alignment/{sample}.bam"
    threads: config["threads"]
    params:
        outbase = "genome/" + basename,
        max_mismatch = config["bowtie"]["max_mismatch"],
        min_ins = config["bowtie"]["min_ins"],
        max_ins = config["bowtie"]["max_ins"]
    log:
       "logs/bowtie/bowtie-align-{sample}.log"
    shell: """
        (bowtie -v {params.max_mismatch} -I {params.min_ins} -X {params.max_ins} --fr --nomaqround --best -S -p {threads} --un alignment/unaligned-{wildcards.sample}.fastq {params.outbase} -1 {input.r1} -2 {input.r2} | samtools view -buh -f 0x2 - | samtools sort -T {wildcards.sample} -@ {threads} -o {output.bam} -) &> {log}
        """

rule samtools_index:
    input:
        "alignment/{sample}.bam"
    output:
        "alignment/{sample}.bam.bai"
    log:
        "logs/samtools_index.log"
    shell: """
        samtools index -b {input} 
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
    shell: """
        samtools view {input} | awk 'BEGIN{{FS=OFS="\t"}} ($2==163) || ($2==99) {{print "+", $3, $4-1, $10}} ($2==83) || ($2==147) {{print "-", $3, $4-1, $10}}' > {output}
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
    shell: """
        samtools sort -n -@ {threads} {input.bam} | bedtools bamtobed -bedpe -i stdin > {output}
        """

rule get_fragsizes:
    input:
        "alignment/fragments/{sample}-fragments.bedpe"
    output:
        temp("alignment/fragments/.{sample}-fragsizes.tsv")
    shell: """
        awk -v sample={wildcards.sample} 'BEGIN{{FS=OFS="\t"}}{{print sample , $6-$2}}' {input} > {output}
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

rule midpoint_coverage:
    input:
        bedpe = "alignment/fragments/{sample}-fragments.bedpe",
        chrsizes = config["genome"]["chrsizes"]
    output:
        "coverage/{sample}-midpoint-counts.bedgraph"
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}} width=$6-$2 {{if(width % 2 != 0){{width -= 1}}; mid=$2+width/2; print $1, mid, mid+1, $7}}' {input.bedpe} | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -g {input.chrsizes} -bga > {output}
        """

rule total_coverage:
    input:
        bam = "alignment/{sample}.bam",
        chrsizes = config["genome"]["chrsizes"]
    output:
        "coverage/{sample}-total-counts.bedgraph"
    shell: """
        bedtools genomecov -ibam {input.bam} -g {input.chrsizes} -bga -pc > {output}
        """

rule seq_depth_norm:
    input:
        midpoint = "coverage/{sample}-midpoint-counts.bedgraph",
        total = "coverage/{sample}-total-counts.bedgraph",
        nfrags = "qual_ctrl/fragment_counts.txt"
    output:
        midpoint = "coverage/{sample}-midpoint-CPM.bedgraph",
        total = "coverage/{sample}-total-CPM.bedgraph"
    shell: """
        {wildcards.sample}norm=$(grep {wildcards.sample} {input.nfrags} | cut -f3 -d ' ')
        awk -v norm=${wildcards.sample}norm 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4*1000000/norm}}' {input.midpoint} > {output.midpoint}
        awk -v norm=${wildcards.sample}norm 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4*1000000/norm}}' {input.total} > {output.total}
        """

rule make_window_file:
    input:
        chrsizes = config["genome"]["chrsizes"]
    output:
        temp("genome/windows.bed")
    params:
        wsize = config["corr-binsize"]
    shell: """
        bedtools makewindows -g {input.chrsizes} -w {params.wsize} | sort -k1,1 -k2,2n > {output}
        """

rule map_to_windows:
    input:
        bed = "genome/windows.bed",
        bedgraph = "coverage/{sample}-midpoint-CPM.bedgraph"
    output:
        temp("coverage/.{sample}-maptowindow.tsv")
    shell: """
        bedtools map -c 4 -o sum -a {input.bed} -b {input.bedgraph} | cut -f4 > {output}
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

                
        






#rule deeptools_multibw_summary:
#    input:
#        expand("coverage/bigwig/{sample}-midpoint-coverage-RPM.bw", sample=SAMPLES)
#    output:
#        matrix = "correlations/midpoint-coverage-RPM.npz",
#        values = "correlations/midpoint-coverage-RPM.tsv"
#    params:
#        labels = list(config["barcodes"].keys()),
#        binsize = config["corr-binsize"]
#    threads: config["threads"]
#    log: "logs/deeptools/multibw_summary.log"
#    shell: """
#        (multiBigwigSummary bins -b {input} -out {output.matrix} --labels {params.labels} -bs {params.binsize} -p {threads} --outRawCounts {output.values}) &> {log}
#        """









