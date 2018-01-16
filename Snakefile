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

localrules: all,
            make_barcode_file,
            bowtie_build,
            samtools_index,
            # build_nucwave_input,
            # get_fragsizes,
            # cat_fragsizes,
            # make_fragments_table,
            # cat_windows,
            # cat_perbase_depth,
            # group_bam_for_danpos,

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
        expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all","passing"]) if sisamples else [],
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-mnase-{{status}}-window-{{windowsize}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]) + expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-mnase-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]) if sisamples else expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-mnase-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]),
        #datavis
        # expand(expand("datavis/{{annotation}}/spikenorm/mnase-{{annotation}}-spikenorm-{{status}}_{condition}-v-{control}-{{readtype}}-{{plot}}-bysample.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), annotation=config["annotations"], readtype=["midpoint","wholefrag"], status=["all","passing"], plot=["heatmap","metagene"]) + expand(expand("datavis/{{annotation}}/libsizenorm/mnase-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}-{{readtype}}-{{plot}}-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], readtype=["midpoint","wholefrag"], status=["all","passing"], plot=["heatmap", "metagene"]) if sisamples else expand(expand("datavis/{{annotation}}/libsizenorm/mnase-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}-{{readtype}}-{{plot}}-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], readtype=["midpoint","wholefrag"], status=["all","passing"], plot=["heatmap","metagene"])

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
    shell: """
        (mkdir -p qual_ctrl/fastqc/raw) &> {log}
        (fastqc -a {input.adapters} --nogroup --extract -t {threads} -o qual_ctrl/fastqc/raw {input.r1} {input.r2}) &>> {log}
        """

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
rule fastqc_processed:
    input:
        "fastq/{fqtype}/{sample}-{fqtype}.{read}.fastq.gz",
    params:
        adapter= lambda wildcards: SAMPLES[wildcards.sample]["barcode"]
    output:
        "qual_ctrl/fastqc/{fqtype}/{sample}-{fqtype}.{read}_fastqc/fastqc_data.txt",
    threads: config["threads"]
    log: "logs/fastqc_processed/fastqc_processed-{sample}-{fqtype}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.fqtype}) &> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/{wildcards.fqtype} {input}) &>> {log}
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
        kmer = 'qual_ctrl/fastqc/kmer_content.tsv'
    output:
        seq_len_dist = 'qual_ctrl/fastqc/sequence_length_distribution.svg',
        per_tile = 'qual_ctrl/fastqc/per_tile_quality.svg',
        per_base_qual = 'qual_ctrl/fastqc/per_base_quality.svg',
        per_base_seq = 'qual_ctrl/fastqc/per_base_sequence_content.svg',
        per_seq_gc = 'qual_ctrl/fastqc/per_sequence_gc.svg',
        per_seq_qual = 'qual_ctrl/fastqc/per_sequence_quality.svg',
        adapter_content = 'qual_ctrl/fastqc/adapter_content.svg',
        seq_dup = 'qual_ctrl/fastqc/sequence_duplication_levels.svg',
        kmer = 'qual_ctrl/fastqc/kmer_content.svg',
    script: "scripts/fastqc_summary.R"

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
    params:
        filterprefix = lambda wildcards: config["combinedgenome"]["spikein_prefix"] if wildcards.species==config["combinedgenome"]["experimental_prefix"] else config["combinedgenome"]["experimental_prefix"],
    output:
        "alignment/{sample}_{species}only.bam"
    threads: config["threads"]
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -h {input.bam} $(grep {wildcards.species} {input.chrsizes} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | grep -v -e 'SN:{params.filterprefix}' | sed 's/{wildcards.species}//g' | samtools view -bh -@ {threads} -o {output} -) &> {log}
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
        "qual_ctrl/{status}/{condition}-v-{control}-mnase-{status}-window-{windowsize}-{norm}-correlations.svg"
    params:
        pcount = 0.1,
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

rule plot_heatmaps:
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

rule plot_metagenes:
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
