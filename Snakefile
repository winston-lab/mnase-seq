#!/usr/bin/env python

configfile: "config.yaml"

SAMPLES = config["barcodes"]

localrules: all,
            make_barcode_file,
            bowtie_build,
            samtools_index,
            build_nucwave_input

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
        "qual_ctrl/frag-size-hist.png",
        #expand("nucwave/{sample}/{sample}.historeadsize.wig", sample=SAMPLES),
        expand("nucwave/{sample}/{sample}_depth_wl_trimmed_PE.wig", sample=SAMPLES),
        #expand("coverage/{sample}-midpoint-coverage-RPM.bedgraph", sample=SAMPLES),
        "correlations/midpoint-coverage-RPM.npz",
        expand("coverage/{sample}-midpoint-coverage-smoothed-RPM.bedgraph", sample=SAMPLES),
        expand("coverage/{sample}-total-coverage-RPM.bedgraph", sample=SAMPLES)

rule make_barcode_file:
    output:
        "barcodes.tsv"
    params:
        bc = config["barcodes"]
    run:
        with open(output[0], "w") as out:
            for x in params.bc:
                out.write(('\t'.join((x, params.bc[x]))+'\n'))

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
# for single end testing:
    #shell: """
    #    (fastq-multx -B {input.barcodes} -b {input.r1} -o fastq/%.r1.fastq.gz) &> {log}
    #    """
# note: we could enforce a minimum quality score for the barcode, but if the quality is that bad at the first bases then the data are probably unusable anyway
# note: for the inline barcodes on both reads, fastq-multx only removes the barcode

#here we use cutadapt to remove the barcode for read2, and cut the A tail and do quality trimming for both reads
#we would like to use --nextseq-trim for 2-color quality trimming, but it doesn't do read2 currently
# note: this assumes ascii(phred+33) quality scores, which I think is standard for modern Illumina sequencers
# THE MINIMUM LENGTH REQUIREMENT IS IMPORTANT TO SANITIZE THE OUTPUT FOR BOWTIE 1
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
# for single end testing:
    #shell: """
    #    (cutadapt -u 1 --nextseq-trim={params.qual_cutoff} -o {output.r1} {input.r1}) &> {log}
    #    """

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

#build bowtie index for given genome
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

#in christine's paper, Burak uses -m 10 --best
rule bowtie:
    input:
        expand("genome/" + basename + ".{n}.ebwt", n=[1,2,3,4]),
        expand("genome/" + basename + ".rev.{n}.ebwt", n=[1,2]),
        r1 = "fastq/trimmed/{sample}-trim.r1.fastq",
        r2 = "fastq/trimmed/{sample}-trim.r2.fastq"
    output:
        bam ="alignment/{sample}.bam"
        #unaligned="alignment/unaligned-{sample}.fastq"
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
#for single end testing
    #shell: """
    #    (bowtie -v {params.max_mismatch} --nomaqround --best -S -p {threads} {params.outbase} {input.r1} | samtools view -buh  - | samtools sort -T {wildcards.sample} -@ {threads} -o {output} -) &> {log}
    #    """

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
        (plotCoverage -n 50000 --outRawCounts {output.table} -p {threads} -v --extendReads --minFragmentLength {params.min_ins} --maxFragmentLength {params.max_ins} -b {input.bam} -o {output.fig}) &> {log}
        """   

rule deeptools_fragsize:
    input:
        bam = expand("alignment/{sample}.bam", sample=SAMPLES),
        index = expand("alignment/{sample}.bam.bai", sample=SAMPLES)
    output:
        hist = "qual_ctrl/frag-size-hist.png"
    log: "logs/deeptools/deeptools_fragsize.log"
    threads: config["threads"]
    shell: """
        (bamPEFragmentSize -b {input.bam} -hist {output.hist} -p {threads}) &> {log} 
        """

rule build_nucwave_input:
    input:
       "alignment/{sample}.bam" 
    output:
       "alignment/{sample}.bowtie"
    shell: """
        samtools view {input} | awk 'BEGIN{{FS=OFS="\t"}} ($2==163) || ($2==99) {{print "+", $3, $4-1, $10}} ($2==83) || ($2==147) {{print "-", $3, $4-1, $10}}' > {output}
        """
#for single end testing:
#    shell: """
#        samtools view {input} | awk 'BEGIN{{FS=OFS="\t"}} ($2==0) {{print "+", $3, $4-1, $10}} ($2==16) {{print "-", $3, $4-1, $10}}' > {output}
#        """

rule nucwave:
    input:
        fasta = config["genome"]["fasta"],
        alignment = "alignment/{sample}.bowtie"
    output:
        "nucwave/{sample}/{sample}_cut_p.wig",
        "nucwave/{sample}/{sample}_cut_m.wig",
#for SE #"nucwave/{sample}/{sample}_depth_p.wig",
        #"nucwave/{sample}/{sample}_depth_m.wig",
        #"nucwave/{sample}/{sample}_depth_p_wl.wig",
        #"nucwave/{sample}/{sample}_depth_m_wl.wig",
        #"nucwave/{sample}/{sample}_depth_c.wig",
        #"nucwave/{sample}/{sample}_depth_c_wl.wig",
        #"nucwave/{sample}/{sample}_depth_c_wl_norm.wig",
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
#for single end testing:
    #shell: """
    #    (python scripts/nucwave_sr.py -w -o nucwave/{wildcards.sample} -g {input.fasta} -a {input.alignment} -p {wildcards.sample}) &> {log}
    #    """

rule midpoint_coverage:
    input:
        bam = "alignment/{sample}.bam",
        index = "alignment/{sample}.bam.bai"
    output:
        "coverage/{sample}-midpoint-coverage-RPM.bedgraph"
    params:
        minsize = config["bowtie"]["min_ins"],
        maxsize = config["bowtie"]["max_ins"]
    threads: config["threads"]
    log: "logs/deeptools/{sample}-midpoint-coverage-RPM.log"
    shell: """
        (bamCoverage -b {input.bam} -o {output} -of bedgraph --MNase -bs 1 -p {threads} --normalizeUsingRPKM --minFragmentLength {params.minsize} --maxFragmentLength {params.maxsize}) &> {log}
        """

rule midpoint_coverage_smoothed:
    input:
        bam = "alignment/{sample}.bam",
        index = "alignment/{sample}.bam.bai"
    output:
        "coverage/{sample}-midpoint-coverage-smoothed-RPM.bedgraph"
    params:
        minsize = config["bowtie"]["min_ins"],
        maxsize = config["bowtie"]["max_ins"],
        smoothwindow = config["coverage"]["smooth_window"]
    threads: config["threads"]
    log: "logs/deeptools/{sample}-midpoint-coverage-smoothed-RPM.bedgraph"
    shell: """
        (bamCoverage -b {input.bam} -o {output} -of bedgraph --MNase -bs 1 -p {threads} --normalizeUsingRPKM --minFragmentLength {params.minsize} --maxFragmentLength {params.maxsize} --smoothLength {params.smoothwindow}) &> {log}
        """

rule total_coverage:
    input:
        bam = "alignment/{sample}.bam",
        index = "alignment/{sample}.bam.bai"
    output:
        "coverage/{sample}-total-coverage-RPM.bedgraph"
    params:
        minsize = config["bowtie"]["min_ins"],
        maxsize = config["bowtie"]["max_ins"]
    threads: config["threads"]
    log: "logs/deeptools/{sample}-total-coverage-RPM.bedgraph"
    shell: """
        (bamCoverage -b {input.bam} -o {output} -of bedgraph --extendReads -bs 1 -p {threads} --normalizeUsingRPKM --minFragmentLength {params.minsize} --maxFragmentLength {params.maxsize}) &> {log}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph = "coverage/{sample}-midpoint-coverage-RPM.bedgraph",
        chrsizes = config["genome"]["chrsizes"]
    output:
        "coverage/bigwig/{sample}-midpoint-coverage-RPM.bw"
    shell: """
        bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}
        """

rule deeptools_multibw_summary:
    input:
        expand("coverage/bigwig/{sample}-midpoint-coverage-RPM.bw", sample=SAMPLES)
    output:
        matrix = "correlations/midpoint-coverage-RPM.npz",
        values = "correlations/midpoint-coverage-RPM.tsv"
    params:
        labels = list(config["barcodes"].keys()),
        binsize = config["corr-binsize"]
    threads: config["threads"]
    log: "logs/deeptools/multibw_summary.log"
    shell: """
        (multiBigwigSummary bins -b {input} -out {output.matrix} --labels {params.labels} -bs {params.binsize} -p {threads} --outRawCounts {output.values}) &> {log}
        """









