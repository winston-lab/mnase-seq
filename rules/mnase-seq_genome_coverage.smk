#!/usr/bin/env python

rule get_fragment_lengths:
    input:
        expand("alignment/{sample}_{species}only.bam", sample=SAMPLES, species=config["combinedgenome"]["experimental_prefix"])
    params:
        header = "\t".join(["fragsize"] + list(SAMPLES.keys()))
    output:
        "qual_ctrl/all/fragment_length_distributions.tsv"
    threads: config["threads"]
    run:
        bam = input[0]
        shell("""samtools view {bam} | cut -f9 | sed 's/-//g' | sort -k1,1n -S 50% --parallel {threads} | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $2, $1}}' > {output}""")
        for bam in input[1:]:
            shell("""join -1 1 -2 2 -t $'\t' -e 0 -a 1 -a 2 {output} <(samtools view {bam} | cut -f9 | sed 's/-//g' | sort -k1,1n -S 50% --parallel {threads} | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}') > qual_ctrl/all/.frag_length.temp; mv qual_ctrl/all/.frag_length.temp {output}""")
        shell("""sed -i "1i {params.header}" {output}""")

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
        bedpe = lambda wc: "alignment/fragments/" + wc.sample + "-" + config["combinedgenome"]["experimental_prefix"] + "fragments.bedpe" if wc.counttype=="counts" else "alignment/fragments/" + wc.sample + "-" + config["combinedgenome"]["spikein_prefix"] + "fragments.bedpe",
        chrsizes = lambda wc: config["genome"]["chrsizes"] if wc.counttype=="counts" else config["genome"]["sichrsizes"]
    params:
        prefix = lambda wc: config["combinedgenome"]["experimental_prefix"] if wc.counttype=="counts" else config["combinedgenome"]["spikein_prefix"]
    output:
        "coverage/{counttype,counts|sicounts}/{sample}-mnase-midpoint-{counttype}.bedgraph"
    log: "logs/midpoint_coverage/midpoint_coverage-{sample}-{counttype}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} {{width=$6-$2}} {{(width % 2 != 0)? (mid=(width+1)/2+$2) : ((rand()<0.5)? (mid=width/2+$2) : (mid=width/2+$2+1))}} {{print $1, mid, mid+1, $7}}' {input.bedpe} | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -g {input.chrsizes} -bga | sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule whole_fragment_coverage:
    input:
        bam = lambda wc: "alignment/" + wc.sample + "_" + config["combinedgenome"]["experimental_prefix"] + "only.bam" if wc.counttype=="counts" else "alignment/" + wc.sample + "_" + config["combinedgenome"]["spikein_prefix"] + "only.bam",
    output:
        "coverage/{counttype}/{sample}-mnase-wholefrag-{counttype}.bedgraph"
    wildcard_constraints:
        counttype="counts|sicounts"
    log : "logs/total_coverage/total_coverage-{sample}-{counttype}.log"
    shell: """
        (bedtools genomecov -ibam {input.bam} -bga -pc | sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule normalize:
    input:
        coverage = "coverage/counts/{sample}-mnase-{readtype}-counts.bedgraph",
        fragcounts = lambda wc: "coverage/counts/" + wc.sample + "-mnase-midpoint-counts.bedgraph" if wc.norm=="libsizenorm" else "coverage/sicounts/" + wc.sample + "-mnase-midpoint-sicounts.bedgraph"
    params:
        scalefactor = lambda wc: config["spikein-pct"] if wc.norm=="spikenorm" else 1
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
        chrsizes = lambda wc: config["genome"]["sichrsizes"] if wc.norm=="sicounts" else config["genome"]["chrsizes"]
    output:
        "coverage/{norm}/{sample}-mnase-{readtype}-{norm}.bw"
    log : "logs/bg_to_bw/bg_to_bw-{sample}-{readtype}-{norm}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
        """

rule smoothed_midpoint_coverage:
    input:
        "coverage/{norm}/{sample}-mnase-midpoint-{norm}.bw"
    output:
        "coverage/{norm}/{sample}-mnase-midpoint_smoothed-{norm}.bw"
    params:
        bandwidth = config["smooth_bandwidth"]
    log: "logs/smoothed_midpoint_coverage/smooth_midpoint_coverage-{sample}-{norm}.log"
    shell: """
        (python scripts/smooth_midpoint_coverage.py -b {params.bandwidth} -i {input} -o {output}) &> {log}
        """

