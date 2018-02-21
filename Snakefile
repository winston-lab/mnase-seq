#!/usr/bin/env python
import itertools

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

FIGURES = config["figures"]

localrules: all,
            make_barcode_file,
            bowtie_build,
            samtools_index,
            cat_matrices,
            group_bam_for_danpos,

rule all:
    input:
        #fastqc
        'qual_ctrl/fastqc/per_base_quality.svg',
        #demultiplex
        expand("fastq/{sample}.{read}.fastq.gz", sample=SAMPLES, read=["r1","r2"]),
        #alignment
        expand("alignment/{sample}.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{counttype}/{sample}-mnase-{readtype}-{counttype}.bedgraph", sample=SAMPLES, readtype=["midpoint","wholefrag"], counttype=COUNTTYPES),
        expand("coverage/{norm}/{sample}-mnase-{readtype}-{norm}.bedgraph", norm=NORMS, sample=SAMPLES, readtype=["midpoint","wholefrag"]),
        expand("coverage/{norm}/{sample}-mnase-midpoint_smoothed-{norm}.bw", norm=NORMS, sample=SAMPLES),
        #quality controls
        "qual_ctrl/read_processing-loss.svg",
        "qual_ctrl/all/fragment_length_distributions.tsv",
        expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all","passing"]) if sisamples else [],
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-mnase-{{status}}-window-{{windowsize}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]) +
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-mnase-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]) if sisamples else expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-mnase-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]),
        #datavis
        # expand("datavis/{figure}/{norm}/allsamples-allannotations-{readtype}-{norm}.tsv.gz", figure=FIGURES, norm=NORMS, readtype=["midpoint", "wholefrag"]),
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]) +
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]) if sisamples else
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]),
        #call nucleosomes 
        expand("nucleosome_calling/{condition}-v-{control}/reference_positions.xls", zip, condition=conditiongroups, control=controlgroups),
        expand("nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_dyad-shift-histogram.svg", zip, condition=conditiongroups, control=controlgroups)

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

#cutadapt doesn't demultiplex paired end, so we use a different script
#NOTE: as of v1.15, cutadapt now demultiplexes paired end, could try it out
# this would allow us to throw out reads where adapter isn't found in both reads
rule demultiplex:
    input:
        r1 = config["fastq"]["r1"],
        r2 = config["fastq"]["r2"],
        barcodes = "fastq/barcodes.tsv"
    output:
        r1 = expand("fastq/{sample}.r1.fastq.gz", sample=["unmatched"] + list(SAMPLES.keys())),
        r2 = expand("fastq/{sample}.r2.fastq.gz", sample=["unmatched"] + list(SAMPLES.keys()))
    log:
        "logs/demultiplex.log"
    shell: """
       (fastq-multx -B {input.barcodes} -b {input.r1} {input.r2} -o fastq/%.r1.fastq.gz -o fastq/%.r2.fastq.gz ) &> {log}
        """
# note: fastq-multx only remove the barcode on read 1. The barcode on read 2 is removed with cutadapt.

#fastQC on raw (demultiplexed) sequencing data
rule fastqc_raw:
    input:
        r1 = "fastq/{sample}.r1.fastq.gz",
        r2 = "fastq/{sample}.r2.fastq.gz",
        adapters = "fastq/barcodes.tsv"
    output:
        "qual_ctrl/fastqc/raw/{sample}.r1_fastqc/fastqc_data.txt",
        "qual_ctrl/fastqc/raw/{sample}.r2_fastqc/fastqc_data.txt",
    threads : config["threads"]
    log : "logs/fastqc_raw/fastqc_raw-{sample}.log"
    run:
        if wildcards.sample=="unmatched":
            shell("""(mkdir -p qual_ctrl/fastqc/raw) &> {log};
                    (fastqc -a {input.adapters} --nogroup --extract -t {threads} -o qual_ctrl/fastqc/raw {input.r1}) &>> {log};
                    (fastqc -a {input.adapters} --nogroup --extract -t {threads} -o qual_ctrl/fastqc/raw {input.r2}) &>> {log}""")
        else:
            adapter = SAMPLES[wildcards.sample]["barcode"]
            shell("""(mkdir -p qual_ctrl/fastqc/raw) &> {log};
                    (fastqc -a <(echo -e "adapter\t{adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/raw {input.r1}) &>> {log};
                    (fastqc -a <(echo -e "adapter\t{adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/raw {input.r2}) &>> {log}""")

# cutadapt:
#    remove barcode from read 2
#    trim the 'A' tail off
#    do quality trimming for both reads
#        - ideally, we would use --nextseq-trim for 2-color quality trimming instead of -q
#            - however, --nextseq-trim currently doesn't trim read 2
#    note: the minimum length requirement (trimmed read >= 5nt) is to sanitize the output for bowtie 1
#    note: the maximum length requirement is to discard reads in which the barcode isn't found in read 2
rule cutadapt:
    input:
        r1 = "fastq/{sample}.r1.fastq.gz",
        r2 = "fastq/{sample}.r2.fastq.gz"
    output:
        r1 = "fastq/cleaned/{sample}-cleaned.r1.fastq.gz",
        r2 = "fastq/cleaned/{sample}-cleaned.r2.fastq.gz",
        log = "logs/cutadapt/cutadapt-{sample}.log"
    params:
        qual_cutoff = config["cutadapt"]["qual_cutoff"],
        adapter = lambda wildcards : SAMPLES[wildcards.sample]["barcode"]+"T",
        max_len = lambda wildcards: config["read-length"] - len(SAMPLES[wildcards.sample]["barcode"]+"T"),
    # threads: config["threads"]
    shell: """
        (cutadapt -e 0.15 -u 1 -G ^{params.adapter} -q {params.qual_cutoff} --minimum-length 5 --maximum-length {params.max_len} -o {output.r1} -p {output.r2} {input.r1} {input.r2}) &> {output.log}
        """

#fastqc for cleaned, aligned, and unaligned reads
#do the two reads sequentially to avoid fastqc bugs and problems with protected files (writing to same directory)
rule fastqc_processed:
    input:
        r1 = "fastq/{fqtype}/{sample}-{fqtype}.r1.fastq.gz",
        r2 = "fastq/{fqtype}/{sample}-{fqtype}.r2.fastq.gz",
    params:
        adapter= lambda wildcards: SAMPLES[wildcards.sample]["barcode"]
    output:
        "qual_ctrl/fastqc/{fqtype}/{sample}-{fqtype}.r1_fastqc/fastqc_data.txt",
        "qual_ctrl/fastqc/{fqtype}/{sample}-{fqtype}.r2_fastqc/fastqc_data.txt",
    threads: config["threads"]
    log: "logs/fastqc_processed/fastqc_processed-{sample}-{fqtype}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.fqtype}) &> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/{wildcards.fqtype} {input.r1}) &>> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/{wildcards.fqtype} {input.r2}) &>> {log}
        """

rule read_processing_numbers:
    input:
        adapter = expand("logs/cutadapt/cutadapt-{sample}.log", sample=SAMPLES),
        align = expand("logs/align/align-{sample}.log", sample=SAMPLES),
    output:
        "qual_ctrl/read_processing_summary.tsv"
    log: "logs/read_processing_summary.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped" > {output}) &> {log}""")
        for sample, adapter, align in zip(SAMPLES.keys(), input.adapter, input.align):
            shell("""(grep -e "Total read pairs processed:" -e "Pairs written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &> {log}""")
            shell("""(grep -e "^Reported" {align} | awk '{{print $2}}' >> {output}) &> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing_summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing-loss.svg",
    script: "scripts/processing_summary.R"

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
rule align:
    input:
        expand(config["bowtie"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".{num}.ebwt", num=[1,2,3,4]) if sisamples else expand(config["bowtie"]["index-path"] + "/" + config["genome"]["name"] + ".{num}.ebwt", num=[1,2,3,4]),
        expand(config["bowtie"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".rev.{num}.ebwt", num=[1,2]) if sisamples else expand(config["bowtie"]["index-path"] + "/" + config["genome"]["name"] + ".rev.{num}.ebwt", num=[1,2]),
        r1 = "fastq/cleaned/{sample}-cleaned.r1.fastq.gz",
        r2 = "fastq/cleaned/{sample}-cleaned.r2.fastq.gz"
    params:
        idx_path = config["bowtie"]["index-path"],
        basename = config["combinedgenome"]["name"] if sisamples else config["genome"]["name"],
        max_mismatch = config["bowtie"]["max_mismatch"],
        min_ins = config["bowtie"]["min_ins"],
        max_ins = config["bowtie"]["max_ins"]
    output:
        bam ="alignment/{sample}.bam",
        aligned_1_gz = "fastq/aligned/{sample}-aligned.r1.fastq.gz",
        aligned_2_gz = "fastq/aligned/{sample}-aligned.r2.fastq.gz",
        unaligned_1_gz = "fastq/unaligned/{sample}-unaligned.r1.fastq.gz",
        unaligned_2_gz = "fastq/unaligned/{sample}-unaligned.r2.fastq.gz",
        log = "logs/align/align-{sample}.log"
    threads: config["threads"]
    shell: """
        (bowtie -v {params.max_mismatch} -I {params.min_ins} -X {params.max_ins} --fr --nomaqround --best -S -p {threads} --al fastq/aligned/{wildcards.sample}-aligned.fastq --un fastq/unaligned/{wildcards.sample}-unaligned.fastq {params.idx_path}/{params.basename} -1 {input.r1} -2 {input.r2} | samtools view -buh -f 0x2 - | samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {output.log}
        (pigz -f fastq/*/{wildcards.sample}-*aligned_*.fastq) &>> {output.log}
        (mv fastq/aligned/{wildcards.sample}-aligned_1.fastq.gz fastq/aligned/{wildcards.sample}-aligned.r1.fastq.gz) &>> {output.log}
        (mv fastq/aligned/{wildcards.sample}-aligned_2.fastq.gz fastq/aligned/{wildcards.sample}-aligned.r2.fastq.gz) &>> {output.log}
        (mv fastq/unaligned/{wildcards.sample}-unaligned_1.fastq.gz fastq/unaligned/{wildcards.sample}-unaligned.r1.fastq.gz) &>> {output.log}
        (mv fastq/unaligned/{wildcards.sample}-unaligned_2.fastq.gz fastq/unaligned/{wildcards.sample}-unaligned.r2.fastq.gz) &>> {output.log}
        """

rule fastqc_aggregate:
    input:
        raw = expand("qual_ctrl/fastqc/raw/{sample}.{read}_fastqc/fastqc_data.txt", sample=["unmatched"] + list(SAMPLES.keys()), read=["r1", "r2"]),
        cleaned = expand("qual_ctrl/fastqc/cleaned/{sample}-cleaned.{read}_fastqc/fastqc_data.txt", sample=SAMPLES, read=["r1","r2"]),
        aligned = expand("qual_ctrl/fastqc/aligned/{sample}-aligned.{read}_fastqc/fastqc_data.txt", sample=SAMPLES, read=["r1","r2"]),
        unaligned = expand("qual_ctrl/fastqc/unaligned/{sample}-unaligned.{read}_fastqc/fastqc_data.txt", sample=SAMPLES, read=["r1","r2"]),
    output:
        'qual_ctrl/fastqc/per_base_quality.tsv',
        'qual_ctrl/fastqc/per_tile_quality.tsv',
        'qual_ctrl/fastqc/per_sequence_quality.tsv',
        'qual_ctrl/fastqc/per_base_sequence_content.tsv',
        'qual_ctrl/fastqc/per_sequence_gc.tsv',
        'qual_ctrl/fastqc/per_base_n.tsv',
        'qual_ctrl/fastqc/sequence_length_distribution.tsv',
        'qual_ctrl/fastqc/sequence_duplication_levels.tsv',
        'qual_ctrl/fastqc/adapter_content.tsv',
        'qual_ctrl/fastqc/kmer_content.tsv'
    run:
        shell("rm -f {output}")
        #for each statistic
        for outpath, stat, header in zip(output, ["Per base sequence quality", "Per tile sequence quality", "Per sequence quality scores", "Per base sequence content", "Per sequence GC content", "Per base N content", "Sequence Length Distribution", "Total Deduplicated Percentage", "Adapter Content", "Kmer Content"], ["base\tmean\tmedian\tlower_quartile\tupper_quartile\tten_pct\tninety_pct\tsample\tstatus", "tile\tbase\tmean\tsample\tstatus",
        "quality\tcount\tsample\tstatus", "base\tg\ta\tt\tc\tsample\tstatus", "gc_content\tcount\tsample\tstatus", "base\tn_count\tsample\tstatus", "length\tcount\tsample\tstatus", "duplication_level\tpct_of_deduplicated\tpct_of_total\tsample\tstatus", "position\tpct\tsample\tstatus",
        "sequence\tcount\tpval\tobs_over_exp_max\tmax_position\tsample\tstatus" ]):
            for input_type in ["raw", "cleaned", "aligned", "unaligned"]:
                sample_id_list = ["_".join(x) for x in itertools.product(["unmatched"]+list(SAMPLES.keys()), ["r1", "r2"])] if input_type=="raw" else ["_".join(x) for x in itertools.product(SAMPLES.keys(), ["r1", "r2"])]
                for sample_id, fqc in zip(sample_id_list, input[input_type]):
                    if sample_id in ["unmatched_r1", "unmatched_r2"] and stat=="Adapter Content":
                        shell("""awk 'BEGIN{{FS=OFS="\t"}} /{stat}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{m=$2;for(i=2;i<=NF-2;i++)if($i>m)m=$i; print $1, m, "{sample_id}", "{input_type}"}}' {fqc} | tail -n +2 >> {outpath}""")
                    else:
                        shell("""awk 'BEGIN{{FS=OFS="\t"}} /{stat}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{print $0, "{sample_id}", "{input_type}"}}' {fqc} | tail -n +2 >> {outpath}""")
            shell("""sed -i "1i {header}" {outpath}""")

rule plot_fastqc_summary:
    input:
        seq_len_dist = 'qual_ctrl/fastqc/sequence_length_distribution.tsv',
        per_tile = 'qual_ctrl/fastqc/per_tile_quality.tsv',
        per_base_qual = 'qual_ctrl/fastqc/per_base_quality.tsv',
        per_base_seq = 'qual_ctrl/fastqc/per_base_sequence_content.tsv',
        per_base_n = 'qual_ctrl/fastqc/per_base_n.tsv',
        per_seq_gc = 'qual_ctrl/fastqc/per_sequence_gc.tsv',
        per_seq_qual = 'qual_ctrl/fastqc/per_sequence_quality.tsv',
        adapter_content = 'qual_ctrl/fastqc/adapter_content.tsv',
        seq_dup = 'qual_ctrl/fastqc/sequence_duplication_levels.tsv',
        # kmer = 'qual_ctrl/fastqc/kmer_content.tsv'
    output:
        seq_len_dist = 'qual_ctrl/fastqc/sequence_length_distribution.svg',
        per_tile = 'qual_ctrl/fastqc/per_tile_quality.svg',
        per_base_qual = 'qual_ctrl/fastqc/per_base_quality.svg',
        per_base_seq = 'qual_ctrl/fastqc/per_base_sequence_content.svg',
        per_seq_gc = 'qual_ctrl/fastqc/per_sequence_gc.svg',
        per_seq_qual = 'qual_ctrl/fastqc/per_sequence_quality.svg',
        adapter_content = 'qual_ctrl/fastqc/adapter_content.svg',
        seq_dup = 'qual_ctrl/fastqc/sequence_duplication_levels.svg',
        # kmer = 'qual_ctrl/fastqc/kmer_content.svg',
    script: "scripts/fastqc_summary.R"

#the index is required to use region arguments in samtools view to separate the species
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
    params:
        filterprefix = lambda wildcards: config["combinedgenome"]["spikein_prefix"] if wildcards.species==config["combinedgenome"]["experimental_prefix"] else config["combinedgenome"]["experimental_prefix"],
    output:
        "alignment/{sample}_{species}only.bam"
    threads: config["threads"]
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -h {input.bam} $(grep {wildcards.species} {input.chrsizes} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | grep -v -e 'SN:{params.filterprefix}' | sed 's/{wildcards.species}//g' | samtools view -bh -@ {threads} -o {output} -) &> {log}
        """

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
        bedpe = lambda wildcards: "alignment/fragments/" + wildcards.sample + "-" + config["combinedgenome"]["experimental_prefix"] + "fragments.bedpe" if wildcards.counttype=="counts" else "alignment/fragments/" + wildcards.sample + "-" + config["combinedgenome"]["spikein_prefix"] + "fragments.bedpe",
        chrsizes = lambda wildcards: config["genome"]["chrsizes"] if wildcards.counttype=="counts" else config["genome"]["sichrsizes"]
    params:
        prefix = lambda wildcards: config["combinedgenome"]["experimental_prefix"] if wildcards.counttype=="counts" else config["combinedgenome"]["spikein_prefix"]
    output:
        "coverage/{counttype,counts|sicounts}/{sample}-mnase-midpoint-{counttype}.bedgraph"
    log: "logs/midpoint_coverage/midpoint_coverage-{sample}-{counttype}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} {{width=$6-$2}} {{(width % 2 != 0)? (mid=(width+1)/2+$2) : ((rand()<0.5)? (mid=width/2+$2) : (mid=width/2+$2+1))}} {{print $1, mid, mid+1, $7}}' {input.bedpe} | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -g {input.chrsizes} -bga | sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule get_si_pct:
    input:
        plmin = expand("coverage/counts/{sample}-mnase-midpoint-counts.bedgraph", sample=sisamples),
        SIplmin = expand("coverage/sicounts/{sample}-mnase-midpoint-sicounts.bedgraph", sample=sisamples)
    params:
        group = [v["group"] for k,v in sisamples.items()]
    output:
        "qual_ctrl/all/spikein-counts.tsv"
    log: "logs/get_si_pct.log"
    run:
        shell("rm -f {output}")
        for name, exp, si, g in zip(sisamples.keys(), input.plmin, input.SIplmin, params.group):
            shell("""(echo -e "{name}\t{g}\t" $(awk 'BEGIN{{FS=OFS="\t"; ex=0; si=0}}{{if(NR==FNR){{si+=$4}} else{{ex+=$4}}}} END{{print ex+si, ex, si}}' {si} {exp}) >> {output}) &> {log}""")

rule plot_si_pct:
    input:
        "qual_ctrl/all/spikein-counts.tsv"
    output:
        plot = "qual_ctrl/{status}/{status}-spikein-plots.svg",
        stats = "qual_ctrl/{status}/{status}-spikein-stats.tsv"
    params:
        samplelist = lambda wildcards : [k for k,v in sisamples.items() if v["spikein"]=="y"] if wildcards.status=="all" else [k for k,v in sipassing.items() if v["spikein"]=="y"],
        conditions = config["comparisons"]["spikenorm"]["conditions"],
        controls = config["comparisons"]["spikenorm"]["controls"],
    script: "scripts/plotsipct.R"

rule whole_fragment_coverage:
    input:
        bam = lambda wildcards: "alignment/" + wildcards.sample + "_" + config["combinedgenome"]["experimental_prefix"] + "only.bam" if wildcards.counttype=="counts" else "alignment/" + wildcards.sample + "_" + config["combinedgenome"]["spikein_prefix"] + "only.bam",
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

rule map_to_windows:
  input:
      bg = "coverage/{norm}/{sample}-mnase-midpoint-{norm}.bedgraph",
      chrsizes = config["genome"]["chrsizes"]
  output:
      temp("coverage/{norm}/{sample}_mnase-midpoint-window-{windowsize}-coverage-{norm}.bedgraph")
  shell: """
    bedtools makewindows -g {input.chrsizes} -w {wildcards.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output}
    """

rule join_window_counts:
    input:
        exp = expand("coverage/{{norm}}/{sample}_mnase-midpoint-window-{{windowsize}}-coverage-{{norm}}.bedgraph", sample=SAMPLES),
    output:
        exp = "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz",
    params:
        names = list(SAMPLES.keys())
    log: "logs/join_window_counts/join_window_counts-{norm}.log"
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output.exp}) &> {log}
        """

rule plotcorrelations:
    input:
        "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz"
    output:
        "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-mnase-{status}-window-{windowsize}-{norm}-correlations.svg"
    params:
        pcount = lambda wildcards: 0.01*int(wildcards.windowsize),
        samplelist = plotcorrsamples
    script:
        "scripts/plotcorr.R"

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

rule smoothed_midpoint_coverage:
    input:
        "coverage/{norm}/{sample}-mnase-midpoint-{norm}.bw"
    output:
        "coverage/{norm}/{sample}-mnase-midpoint_smoothed-{norm}.bw"
    params:
        bandwidth = config["smooth_bandwidth"]
    log: "logs/smoothed_midpoint_coverage/smooth_midpoint_coverage-{sample}.log"
    shell: """
        (python scripts/smooth_midpoint_coverage.py -b {params.bandwidth} -i {input} -o {output}) &> {log}
        """

rule compute_matrix:
    input:
        annotation = lambda wildcards: FIGURES[wildcards.figure]["annotations"][wildcards.annotation]["path"],
        bw = "coverage/{norm}/{sample}-mnase-{readtype}-{norm}.bw"
    output:
        dtfile = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{readtype}-{norm}.mat.gz"),
        matrix = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{readtype}-{norm}.tsv"),
        melted = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{readtype}-{norm}-melted.tsv.gz")
    params:
        group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
        refpoint = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["refpoint"],
        upstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["upstream"] + FIGURES[wildcards.figure]["parameters"]["binsize"],
        dnstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["dnstream"] + FIGURES[wildcards.figure]["parameters"]["binsize"],
        binsize = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["binsize"],
        binstat = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["binstat"],
        nan_afterend = lambda wildcards: "--nanAfterEnd" if FIGURES[wildcards.figure]["parameters"]["nan_afterend"] else [],
        anno_label = lambda wildcards: FIGURES[wildcards.figure]["annotations"][wildcards.annotation]["label"]
    threads : config["threads"]
    log: "logs/compute_matrix/compute_matrix-{figure}_{annotation}_{sample}_{readtype}_{norm}.log"
    run:
        shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} {params.nan_afterend} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {params.refpoint} -g {params.group} -s {wildcards.sample} -a {params.anno_label} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wildcards: expand("datavis/{figure}/{norm}/{annotation}_{sample}-{readtype}-{norm}-melted.tsv.gz", annotation=[k for k,v in FIGURES[wildcards.figure]["annotations"].items()], sample=SAMPLES, figure=wildcards.figure, norm=wildcards.norm, readtype=wildcards.readtype)
    output:
        "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-{readtype}-{norm}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{figure}-{readtype}-{norm}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_figures:
    input:
        matrix = "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-{readtype}-{norm}.tsv.gz",
        annotations = lambda wildcards: [v["path"] for k,v in FIGURES[wildcards.figure]["annotations"].items()]
    output:
        heatmap_sample = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bysample.svg",
        heatmap_group = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bygroup.svg",
        meta_sample = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bysample.svg",
        meta_sample_overlay = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bysample-overlay.svg",
        meta_group = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bygroup.svg",
        meta_sample_clust = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bycluster-sample.svg",
        meta_group_clust = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/mnase-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bycluster-group.svg",
    params:
        # abusing snakemake a bit here...using params as output paths since in order to use lambda functions
        annotations_out = lambda wildcards: ["datavis/" + wildcards.figure + "/" + wildcards.norm + "/" + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.status + "/" + wildcards.readtype + "/" + annotation + "_cluster-" + str(cluster) + ".bed" for annotation in FIGURES[wildcards.figure]["annotations"] for cluster in range(1, FIGURES[wildcards.figure]["annotations"][annotation]["n_clusters"]+1)],
        clusters_out = lambda wildcards: ["datavis/" + wildcards.figure + "/" + wildcards.norm + "/" + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.status + "/" + wildcards.readtype + "/" + annotation + ".pdf" for annotation in FIGURES[wildcards.figure]["annotations"]],
        samplelist = plotcorrsamples,
        plottype = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["type"],
        readtype = lambda wildcards: "dyad signal" if wildcards.readtype=="midpoint" else "protection",
        upstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["upstream"],
        dnstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["dnstream"],
        scaled_length = lambda wildcards: 0 if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute" else FIGURES[wildcards.figure]["parameters"]["scaled_length"],
        pct_cutoff = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["pct_cutoff"],
        trim_pct = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["trim_pct"],
        refpointlabel = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["refpointlabel"],
        endlabel = lambda wildcards:  "HAIL SATAN" if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute" else FIGURES[wildcards.figure]["parameters"]["endlabel"],
        cmap = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["heatmap_colormap"],
        sortmethod = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["arrange"],
        cluster_groups = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["cluster_conditions"],
        cluster_five = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["cluster_five"],
        cluster_three = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["cluster_three"],
        k = lambda wildcards: [v["n_clusters"] for k,v in FIGURES[wildcards.figure]["annotations"].items()],
    script:
        "scripts/plot_mnase_figures.R"

rule group_bam_for_danpos:
    input:
        "alignment/{sample}_" + config["combinedgenome"]["experimental_prefix"] + "only.bam"
    output:
        "nucleosome_calling/data/{group}/{sample}.bam"
    shell: """
        cp {input} {output}
        """

#TODO: spike-in normalization with --count
rule danpos:
    input:
        ["nucleosome_calling/data/" + SAMPLES[x]['group'] + "/" + x + ".bam" for x in SAMPLES]
    output:
        "nucleosome_calling/{condition}-v-{control}/nucleosome_calling_data_{condition}-nucleosome_calling_data_{control}.positions.integrative.xls",
        "nucleosome_calling/{condition}-v-{control}/reference_positions.xls",
        "nucleosome_calling/{condition}-v-{control}/diff/nucleosome_calling_data_{condition}-nucleosome_calling_data_{control}.pois_diff.wig",
        "nucleosome_calling/{condition}-v-{control}/pooled/nucleosome_calling_data_{condition}.Fnor.smooth.positions.ref_adjust.xls",
        "nucleosome_calling/{condition}-v-{control}/pooled/nucleosome_calling_data_{condition}.Fnor.smooth.positions.xls",
        "nucleosome_calling/{condition}-v-{control}/pooled/nucleosome_calling_data_{condition}.Fnor.smooth.wig",
        "nucleosome_calling/{condition}-v-{control}/pooled/nucleosome_calling_data_{control}.Fnor.smooth.positions.ref_adjust.xls",
        "nucleosome_calling/{condition}-v-{control}/pooled/nucleosome_calling_data_{control}.Fnor.smooth.positions.xls",
        "nucleosome_calling/{condition}-v-{control}/pooled/nucleosome_calling_data_{control}.Fnor.smooth.wig"
    conda:
        "envs/danpos.yaml"
    log: "logs/danpos/danpos-{condition}-v-{control}.log"
    shell: """
        (python scripts/danpos-2.2.2/danpos.py dpos nucleosome_calling/data/{wildcards.condition}/:nucleosome_calling/data/{wildcards.control}/ -m 1 -e 1 -a 1 -o nucleosome_calling/{wildcards.condition}-v-{wildcards.control}) &> {log}
        """

rule plot_danpos_results:
    input:
        results = "nucleosome_calling/{condition}-v-{control}/nucleosome_calling_data_{condition}-nucleosome_calling_data_{control}.positions.integrative.xls",
    output:
        shift_hist = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_dyad-shift-histogram.svg",
        occupancy_volcano = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_occupancy_volcano.svg",
        fuzziness_volcano = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_fuzziness_volcano.svg",
        occupancy_v_fuzziness = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_occupancy-v-fuzziness-scatter.svg",
        occupancy_violin = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_occupancy-violin.svg",
        occupancy_freqpoly = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_occupancy-freqpoly.svg",
        occupancy_ecdf = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_occupancy-ecdf.svg",
        fuzziness_violin = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_fuzziness-violin.svg",
        fuzziness_freqpoly = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_fuzziness-freqpoly.svg",
        fuzziness_ecdf = "nucleosome_calling/{condition}-v-{control}/{condition}-v-{control}_fuzziness-ecdf.svg"
    script:
        "scripts/nucleosomes.R"
