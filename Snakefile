#!/usr/bin/env python

import os
import re
import itertools
from math import log2

configfile: "config.yaml"

subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

configfile: build_annotations("config.yaml")

SAMPLES = config["samples"]
SISAMPLES = {k:v for k,v in SAMPLES.items() if v["spikein"]}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}
SIPASSING = {k:v for k,v in PASSING.items() if v["spikein"]}

controlgroups = list(itertools.chain(*[d.values() for d in config["comparisons"]["libsizenorm"]]))
conditiongroups = list(itertools.chain(*[d.keys() for d in config["comparisons"]["libsizenorm"]]))

comparisons_si = config["comparisons"]["spikenorm"]
if comparisons_si:
    controlgroups_si = list(itertools.chain(*[d.values() for d in config["comparisons"]["spikenorm"]]))
    conditiongroups_si = list(itertools.chain(*[d.keys() for d in config["comparisons"]["spikenorm"]]))

FIGURES = config["figures"]
QUANT = config["quantification"]

wildcard_constraints:
    sample = "|".join(re.escape(x) for x in list(SAMPLES.keys()) + ["unmatched"]),
    group = "|".join(set(re.escape(v["group"]) for k,v in SAMPLES.items())),
    control = "|".join(set(re.escape(x) for x in controlgroups + (conditiongroups_si if comparisons_si else []) + ["all"])),
    condition = "|".join(set(re.escape(x) for x in conditiongroups + (controlgroups_si if comparisons_si else []) + ["all"])),
    species = "experimental|spikein",
    read_status = "raw|cleaned|aligned|unaligned",
    figure = "|".join(re.escape(x) for x in list(FIGURES.keys()) + list(QUANT.keys())),
    annotation = "|".join(re.escape(x) for x in set(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES] + [QUANT[figure]["annotations"].keys() for figure in QUANT]))),
    status = "all|passing",
    counttype= "counts|sicounts",
    norm = "counts|sicounts|libsizenorm|spikenorm",
    readtype = "midpoint|wholefrag",
    windowsize = "\d+",
    direction = "all|up|unchanged|down",

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

def statuscheck(dict1, dict2):
    return(["passing"] if dict1 == dict2 else ["all", "passing"])

def conditioncheck(conditionlist):
    return(conditionlist if len(conditionlist)==1 else conditionlist + ["all"])

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
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
        expand("qual_ctrl/spikein/mnase-seq_spikein-plots-{status}.svg", status=statuscheck(SISAMPLES, SIPASSING)) if SISAMPLES else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_mnase-seq-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), status=statuscheck(SISAMPLES, SIPASSING), windowsize=config["scatterplot_binsizes"]) if SISAMPLES and comparisons_si else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_mnase-seq-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), status=statuscheck(SAMPLES, SISAMPLES), windowsize=config["scatterplot_binsizes"]),
        #datavis
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-seq_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), figure=FIGURES, readtype=["midpoint","wholefrag"], status=statuscheck(SISAMPLES, SIPASSING)) if config["plot_figures"] and SISAMPLES and comparisons_si else [],
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-seq_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), figure=FIGURES, readtype=["midpoint","wholefrag"], status=statuscheck(SAMPLES, PASSING)) if config["plot_figures"] else [],
        #call nucleosomes
        expand("nucleosome_quantification/{condition}-v-{control}/spikenorm/reference_positions.xls", zip, condition=conditiongroups_si, control=controlgroups_si) if SIPASSING and comparisons_si else [],
        expand("nucleosome_quantification/{condition}-v-{control}/libsizenorm/reference_positions.xls", zip, condition=conditiongroups, control=controlgroups),
        expand("nucleosome_quantification/{condition}-v-{control}/spikenorm/{condition}-v-{control}_spikenorm-dyad-shift-histogram.svg", zip, condition=conditiongroups_si, control=controlgroups_si) if SIPASSING and comparisons_si else [],
        expand("nucleosome_quantification/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_libsizenorm-dyad-shift-histogram.svg", zip, condition=conditiongroups, control=controlgroups),
        #danpos over annotations
        expand(expand("nucleosome_quantification/regions/{{figure}}/spikenorm/{condition}-v-{control}/{{figure}}_{condition}-v-{control}_spikenorm-individual-occupancy-heatmaps.svg", zip, condition=conditiongroups_si, control=controlgroups_si), figure=QUANT) if SIPASSING and comparisons_si else [],
        expand(expand("nucleosome_quantification/regions/{{figure}}/libsizenorm/{condition}-v-{control}/{{figure}}_{condition}-v-{control}_libsizenorm-individual-occupancy-heatmaps.svg", zip, condition=conditiongroups, control=controlgroups), figure=QUANT),
        #differential nucleosome levels over transcripts
        expand("diff_levels/{condition}-v-{control}/spikenorm/{condition}-v-{control}-mnase-seq-results-spikenorm-all.tsv", zip, condition=conditiongroups_si, control=controlgroups_si) if SIPASSING and comparisons_si else [],
        expand("diff_levels/{condition}-v-{control}/libsizenorm/{condition}-v-{control}-mnase-seq-results-libsizenorm-all.tsv", zip, condition=conditiongroups, control=controlgroups),

