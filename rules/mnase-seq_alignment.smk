#!/usr/bin/env python

localrules:
    build_combined_genome,
    bowtie_build,
    index_bam,

basename = "{exp_name}_{exp_fasta}_{si_name}_{si_fasta}".format(exp_name = config["genome"]["name"],
                                                                exp_fasta = os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0],
                                                                si_name = config["spike_in"]["name"],
                                                                si_fasta = os.path.splitext(os.path.basename(config["spike_in"]["fasta"]))[0]) if SISAMPLES else os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0]

rule build_combined_genome:
    input:
        experimental = os.path.abspath(build_annotations(config["genome"]["fasta"])),
        spikein = config["spike_in"]["fasta"] if SISAMPLES else []
    output:
        "{directory}/{bn}.fa".format(directory = os.path.split(os.path.abspath(build_annotations(config["genome"]["fasta"])))[0], bn=basename),
    params:
        exp_name = config["genome"]["name"],
        si_name = config["spike_in"]["name"] if SISAMPLES else []
    log: "logs/build_combined_genome.log"
    shell: """
        (sed 's/>/>{params.exp_name}_/g' {input.experimental} | \
        cat - <(sed 's/>/>{params.si_name}_/g' {input.spikein}) > {output}) &> {log}
        """

rule bowtie_build:
    input:
        "{directory}/{bn}.fa".format(directory = os.path.split(os.path.abspath(build_annotations(config["genome"]["fasta"])))[0], bn=basename) if SISAMPLES else os.path.abspath(build_annotations(config["genome"]["fasta"])),
    output:
        expand(config["bowtie"]["index-path"] + "/{{basename}}.{num}.ebwt", num=[1,2,3,4]),
        expand(config["bowtie"]["index-path"] + "/{{basename}}.rev.{num}.ebwt", num=[1,2]),
    params:
        idx_path = config["bowtie"]["index-path"],
    conda: "../envs/bowtie.yaml"
    log: "logs/bowtie-build_{basename}.log"
    shell: """
        (bowtie-build {input} {params.idx_path}/{wildcards.basename}) &> {log}
        """

#align with Bowtie 1
#in Christine's paper, Burak uses -m 10 --best
rule align:
    input:
        expand("{directory}/{bn}.{num}.ebwt", directory = config["bowtie"]["index-path"],
                                             bn = basename,
                                             num = [1,2,3,4]),
        expand("{directory}/{bn}.rev.{num}.ebwt", directory = config["bowtie"]["index-path"],
                                             bn = basename,
                                             num = [1,2]),
        r1 = "fastq/cleaned/{sample}-cleaned.r1.fastq.gz",
        r2 = "fastq/cleaned/{sample}-cleaned.r2.fastq.gz"
    output:
        bam ="alignment/{sample}_mnase-seq.bam",
        aligned_1_gz = "fastq/aligned/{sample}-aligned.r1.fastq.gz",
        aligned_2_gz = "fastq/aligned/{sample}-aligned.r2.fastq.gz",
        unaligned_1_gz = "fastq/unaligned/{sample}-unaligned.r1.fastq.gz",
        unaligned_2_gz = "fastq/unaligned/{sample}-unaligned.r2.fastq.gz",
        log = "logs/align/align-{sample}.log"
    params:
        idx_path = config["bowtie"]["index-path"],
        max_mismatch = config["bowtie"]["max_mismatch"],
        min_ins = config["bowtie"]["min_ins"],
        max_ins = config["bowtie"]["max_ins"]
    conda: "../envs/bowtie.yaml"
    threads: config["threads"]
    shell: """
        (bowtie -v {params.max_mismatch} -I {params.min_ins} -X {params.max_ins} --fr --nomaqround --best -S -p {threads} --al fastq/aligned/{wildcards.sample}-aligned.fastq --un fastq/unaligned/{wildcards.sample}-unaligned.fastq {params.idx_path}/{basename} -1 {input.r1} -2 {input.r2} | samtools view -buh -f 0x2 - | samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {output.log}
        (pigz -f fastq/*/{wildcards.sample}-*aligned_*.fastq) &>> {output.log}
        (mv fastq/aligned/{wildcards.sample}-aligned_1.fastq.gz fastq/aligned/{wildcards.sample}-aligned.r1.fastq.gz) &>> {output.log}
        (mv fastq/aligned/{wildcards.sample}-aligned_2.fastq.gz fastq/aligned/{wildcards.sample}-aligned.r2.fastq.gz) &>> {output.log}
        (mv fastq/unaligned/{wildcards.sample}-unaligned_1.fastq.gz fastq/unaligned/{wildcards.sample}-unaligned.r1.fastq.gz) &>> {output.log}
        (mv fastq/unaligned/{wildcards.sample}-unaligned_2.fastq.gz fastq/unaligned/{wildcards.sample}-unaligned.r2.fastq.gz) &>> {output.log}
        """

#the index is required to use region arguments in samtools view to separate the species
rule index_bam:
    input:
        "alignment/{sample}_mnase-seq.bam",
    output:
        "alignment/{sample}_mnase-seq.bam.bai",
    log:
        "logs/samtools_index/samtools_index-{sample}.log"
    shell: """
        (samtools index -b {input}) &> {log}
        """

rule bam_separate_species:
    input:
        bam = "alignment/{sample}_mnase-seq.bam",
        bai = "alignment/{sample}_mnase-seq.bam.bai",
        fasta = "{directory}/{bn}.fa".format(directory = os.path.split(build_annotations(config["genome"]["fasta"]))[0], bn=basename) if SISAMPLES else [],
    output:
        "alignment/{sample}_mnase-seq-{species}.bam"
    params:
        keep_prefix = lambda wc: config["genome"]["name"] if wc.species=="experimental" else config["spike_in"]["name"],
        discard_prefix = lambda wc: config["spike_in"]["name"] if wc.species=="experimental" else config["genome"]["name"],
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    threads: config["threads"]
    shell: """
        (samtools view -h {input.bam} $(grep {params.keep_prefix}_ <(faidx {input.fasta} -i chromsizes) | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | grep -v -e 'SN:{params.discard_prefix}_' | sed 's/{params.keep_prefix}_//g' | samtools view -bh -@ {threads} -o {output} -) &> {log}
        """

