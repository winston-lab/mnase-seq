#!/usr/bin/env python
import os
import re
import itertools
from math import log2

configfile: "config.yaml"

SAMPLES = config["samples"]
SISAMPLES = {k:v for k,v in SAMPLES.items() if v["spikein"]}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}
SIPASSING = {k:v for k,v in PASSING.items() if v["spikein"]}

controlgroups = [v for k,v in config["comparisons"]["libsizenorm"].items()]
conditiongroups = [k for k,v in config["comparisons"]["libsizenorm"].items()]

if SIPASSING:
    import csv # for danpos2 normalization function
    controlgroups_si = [v for k,v in config["comparisons"]["spikenorm"].items()]
    conditiongroups_si = [k for k,v in config["comparisons"]["spikenorm"].items()]

FIGURES = config["figures"]
QUANT = config["quantification"]

status_norm_sample_dict = {
    "all":
        {   "libsizenorm" : SAMPLES,
            "spikenorm" : SISAMPLES
        },
    "passing":
        {   "libsizenorm" : PASSING,
            "spikenorm" : SIPASSING
        }
    }

def get_samples(status, norm, groups):
    if "all" in groups:
        return(list(status_norm_sample_dict[status][norm].keys()))
    else:
        return([k for k,v in status_norm_sample_dict[status][norm].items() if v["group"] in groups])

include: "rules/mnase-seq_clean_reads.smk"
include: "rules/mnase-seq_alignment.smk"
include: "rules/mnase-seq_fastqc.smk"
include: "rules/mnase-seq_library_processing_summary.smk"
include: "rules/mnase-seq_genome_coverage.smk"
include: "rules/mnase-seq_sample_similarity.smk"
include: "rules/mnase-seq_datavis.smk"
include: "rules/mnase-seq_quantification.smk"
include: "rules/mnase-seq_differential_occupancy.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules: all

rule all:
    input:
        #fastqc
        'qual_ctrl/fastqc/mnase-seq-per_base_quality.svg',
        #alignment
        expand("alignment/{sample}_mnase-seq.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}_mnase-{readtype}-{norm}.bw", sample=SAMPLES, norm=["counts", "libsizenorm"], readtype=["midpoint","wholefrag"]),
        expand("coverage/{norm}/{sample}_mnase-{readtype}-{norm}.bw", sample=SISAMPLES, norm=["sicounts", "spikenorm"], readtype=["midpoint","wholefrag"]),
        expand("coverage/libsizenorm/{sample}_mnase-midpoint_smoothed-libsizenorm.bw", sample=SAMPLES),
        expand("coverage/spikenorm/{sample}_mnase-midpoint_smoothed-spikenorm.bw", sample=SISAMPLES),
        #quality controls
        "qual_ctrl/read_processing/mnase-seq_read_processing-loss.svg",
        "qual_ctrl/fragment_length_distributions/mnase-seq_fragment_length_distributions.svg",
        expand("qual_ctrl/spikein/mnase-seq_spikein-plots-{status}.svg", status=["all","passing"]) if SISAMPLES else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_mnase-seq-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]) if SISAMPLES else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_mnase-seq-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]),
        #datavis
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-seq_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]) if config["plot_figures"] and SISAMPLES else [],
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-seq_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]) if config["plot_figures"] else [],
        #call nucleosomes
        expand("nucleosome_quantification/{condition}-v-{control}/spikenorm/reference_positions.xls", zip, condition=conditiongroups_si, control=controlgroups_si) if SIPASSING else [],
        expand("nucleosome_quantification/{condition}-v-{control}/libsizenorm/reference_positions.xls", zip, condition=conditiongroups, control=controlgroups),
        expand("nucleosome_quantification/{condition}-v-{control}/spikenorm/{condition}-v-{control}_spikenorm-dyad-shift-histogram.svg", zip, condition=conditiongroups_si, control=controlgroups_si) if SIPASSING else [],
        expand("nucleosome_quantification/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_libsizenorm-dyad-shift-histogram.svg", zip, condition=conditiongroups, control=controlgroups),
        #danpos over annotations
        expand(expand("nucleosome_quantification/regions/{{figure}}/spikenorm/{condition}-v-{control}/{{figure}}_{condition}-v-{control}_spikenorm-individual-occupancy-heatmaps.svg", zip, condition=conditiongroups_si, control=controlgroups_si), figure=QUANT) if SIPASSING else [],
        expand(expand("nucleosome_quantification/regions/{{figure}}/libsizenorm/{condition}-v-{control}/{{figure}}_{condition}-v-{control}_libsizenorm-individual-occupancy-heatmaps.svg", zip, condition=conditiongroups, control=controlgroups), figure=QUANT),
        #differential nucleosome levels over transcripts
        expand("diff_levels/{condition}-v-{control}/spikenorm/{condition}-v-{control}-mnase-seq-results-spikenorm-all.tsv", zip, condition=conditiongroups_si, control=controlgroups_si) if SIPASSING else [],
        expand("diff_levels/{condition}-v-{control}/libsizenorm/{condition}-v-{control}-mnase-seq-results-libsizenorm-all.tsv", zip, condition=conditiongroups, control=controlgroups),

