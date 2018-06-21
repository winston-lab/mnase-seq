#!/usr/bin/env python

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

#fastqc for cleaned, aligned, and unaligned reads
#do the two reads sequentially to avoid fastqc bugs and problems with protected files (writing to same directory)
rule fastqc_processed:
    input:
        r1 = "fastq/{fqtype}/{sample}-{fqtype}.r1.fastq.gz",
        r2 = "fastq/{fqtype}/{sample}-{fqtype}.r2.fastq.gz",
    params:
        adapter= lambda wc: SAMPLES[wc.sample]["barcode"]
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

