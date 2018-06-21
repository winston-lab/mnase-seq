#!/usr/bin/env python

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
        samplelist = lambda wc : sisamples if wc.status=="all" else sipassing,
        conditions = conditiongroups_si if sisamples else [],
        controls = controlgroups_si if sisamples else [],
    script: "scripts/plotsipct.R"

