#!/usr/bin/env python

localrules:
    normalize_genome_coverage,
    bedgraph_to_bigwig,
    smoothed_midpoint_coverage


#bam must be sorted by name for bedpe. We don't do this in the bowtie step since samtools indexing required position-sorted bam.
rule get_fragments:
    input:
        bam = "alignment/{sample}_mnase-seq-{species}.bam" if SISAMPLES else "alignment/{sample}_mnase-seq.bam"
    output:
        "alignment/fragments/{sample}_{species}-fragments.bedpe"
    threads: config["threads"]
    log : "logs/get_fragments/get_fragments_{sample}-{species}.log"
    shell: """
        (samtools sort -n -T .{wildcards.sample}_{wildcards.species} -@ {threads} {input.bam} | bedtools bamtobed -bedpe -i stdin > {output}) &> {log}
        """

rule midpoint_coverage:
    input:
        bedpe = lambda wc: f"alignment/fragments/{wc.sample}_experimental-fragments.bedpe" if wc.counttype=="counts" else f"alignment/fragments/{wc.sample}_spikein-fragments.bedpe",
        fasta = lambda wc: os.path.abspath(build_annotations(config["genome"]["fasta"])) if wc.counttype=="counts" else config["spike_in"]["fasta"]
    output:
        "coverage/{counttype,counts|sicounts}/{sample}_mnase-midpoint-{counttype}.bedgraph"
    log: "logs/midpoint_coverage/midpoint_coverage_{sample}-{counttype}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} {{width=$6-$2}} {{(width % 2 != 0)? (mid=(width+1)/2+$2) : ((rand()<0.5)? (mid=width/2+$2) : (mid=width/2+$2+1))}} {{print $1, mid, mid+1, $7}}' {input.bedpe} | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -g <(faidx {input.fasta} -i chromsizes) -bga | sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule whole_fragment_coverage:
    input:
        bam = lambda wc: f"alignment/{wc.sample}_mnase-seq-experimental.bam" if wc.counttype=="counts" and SISAMPLES else f"alignment/{wc.sample}_mnase-seq.bam" if wc.counttype=="counts" else f"alignment/{wc.sample}_mnase-seq-spikein.bam",
    output:
        "coverage/{counttype}/{sample}_mnase-wholefrag-{counttype}.bedgraph"
    wildcard_constraints:
        counttype="counts|sicounts"
    log : "logs/whole_fragment_coverage/whole_fragment_coverage_{sample}-{counttype}.log"
    shell: """
        (bedtools genomecov -ibam {input.bam} -bga -pc | sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_mnase-{readtype}-counts.bedgraph",
        bam = lambda wc: f"alignment/{wc.sample}_mnase-seq-experimental.bam" if wc.norm=="libsizenorm" and SISAMPLES else f"alignment/{wc.sample}_mnase-seq.bam" if wc.norm=="libsizenorm" else f"alignment/{wc.sample}_mnase-seq-spikein.bam"
    output:
        normalized = "coverage/{norm}/{sample}_mnase-{readtype}-{norm}.bedgraph"
    params:
        scale_factor = lambda wc: config["spike_in"]["proportion"] if wc.norm=="spikenorm" else 1
    wildcard_constraints:
        norm="libsizenorm|spikenorm"
    log: "logs/normalize_genome_coverage/normalize_genome_coverage_{sample}-{norm}-{readtype}.log"
    shell: """
        (awk -v norm_factor=$(samtools view -c {input.bam} | paste -d "" - <(echo "/({params.scale_factor}*1000000)") | bc -l) 'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph = "coverage/{norm}/{sample}_mnase-{readtype}-{norm}.bedgraph",
        fasta = lambda wc: config["spike_in"]["fasta"] if wc.norm=="sicounts" else os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        "coverage/{norm}/{sample}_mnase-{readtype}-{norm}.bw"
    log : "logs/bedgraph_to_bigwig/bedgraph_to_bigwig_{sample}-{readtype}-{norm}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} <(faidx {input.fasta} -i chromsizes) {output}) &> {log}
        """

rule smoothed_midpoint_coverage:
    input:
        "coverage/{norm}/{sample}_mnase-midpoint-{norm}.bw"
    output:
        "coverage/{norm}/{sample}_mnase-midpoint_smoothed-{norm}.bw"
    params:
        bandwidth = config["smooth_bandwidth"]
    conda: "../envs/smooth_coverage.yaml"
    log: "logs/smoothed_midpoint_coverage/smoothed_midpoint_coverage_{sample}-{norm}.log"
    shell: """
        (python scripts/smooth_midpoint_coverage.py -b {params.bandwidth} -i {input} -o {output}) &> {log}
        """

