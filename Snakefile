#!/usr/bin/env python
from math import log2
import itertools

configfile: "config.yaml"

SAMPLES = config["samples"]
sisamples = {k:v for k,v in SAMPLES.items() if v["spikein"] == "y"}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"] == "pass"}
sipassing = {k:v for k,v in PASSING.items() if v["spikein"] == "y"}

# if there are samples with spikein, import csv for danpos2 step
if sipassing:
    import csv

controlgroups = [v for k,v in config["comparisons"]["libsizenorm"].items()]
conditiongroups = [k for k,v in config["comparisons"]["libsizenorm"].items()]

if sipassing:
    controlgroups_si = [v for k,v in config["comparisons"]["spikenorm"].items()]
    conditiongroups_si = [k for k,v in config["comparisons"]["spikenorm"].items()]

COUNTTYPES = ["counts", "sicounts"] if sisamples else ["counts"]
NORMS = ["libsizenorm", "spikenorm"] if sisamples else ["libsizenorm"]

FIGURES = config["figures"]
QUANT = config["quantification"]

localrules: all,
            make_barcode_file,
            bowtie_build,
            samtools_index,
            cat_matrices,
            group_bam_for_danpos,
            danpos_over_annotations,
            cat_danpos_annotations,
            map_counts_to_transcripts, get_transcript_counts,

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

rule all:
    input:
        ##fastqc
        #'qual_ctrl/fastqc/per_base_quality.svg',
        ##alignment
        #expand("alignment/{sample}.bam", sample=SAMPLES),
        ##coverage
        #expand("coverage/{counttype}/{sample}-mnase-{readtype}-{counttype}.bedgraph", sample=SAMPLES, readtype=["midpoint","wholefrag"], counttype=COUNTTYPES),
        #expand("coverage/{norm}/{sample}-mnase-{readtype}-{norm}.bedgraph", norm=NORMS, sample=SAMPLES, readtype=["midpoint","wholefrag"]),
        #expand("coverage/{norm}/{sample}-mnase-midpoint_smoothed-{norm}.bw", norm=NORMS, sample=SAMPLES),
        ##quality controls
        #"qual_ctrl/read_processing-loss.svg",
        #"qual_ctrl/all/fragment_length_distributions.tsv",
        #expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all","passing"]) if sisamples else [],
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-mnase-{{status}}-window-{{windowsize}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]) +
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-mnase-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]) if sisamples else expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-mnase-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all","passing"], windowsize=config["corr-windowsizes"]),
        ##datavis
        #expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]) +
        #expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]) if sisamples else
        #expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]),
        ##call nucleosomes
        #expand("nucleosome_calling/{condition}-v-{control}/spikenorm/reference_positions.xls", zip, condition=conditiongroups_si, control=controlgroups_si) + expand("nucleosome_calling/{condition}-v-{control}/libsizenorm/reference_positions.xls", zip, condition=conditiongroups, control=controlgroups) if sipassing else expand("nucleosome_calling/{condition}-v-{control}/libsizenorm/reference_positions.xls", zip, condition=conditiongroups, control=controlgroups),
        #expand("nucleosome_calling/{condition}-v-{control}/spikenorm/{condition}-v-{control}_spikenorm-dyad-shift-histogram.svg", zip, condition=conditiongroups_si, control=controlgroups_si) + expand("nucleosome_calling/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_libsizenorm-dyad-shift-histogram.svg", zip, condition=conditiongroups, control=controlgroups) if sipassing else expand("nucleosome_calling/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_libsizenorm-dyad-shift-histogram.svg", zip, condition=conditiongroups, control=controlgroups),
        ##danpos over annotations
        #expand(expand("nucleosome_calling/regions/{{figure}}/spikenorm/{condition}-v-{control}/{{figure}}_{condition}-v-{control}_spikenorm-individual-occupancy-heatmaps.svg", zip, condition=conditiongroups_si, control=controlgroups_si), figure=QUANT) + expand(expand("nucleosome_calling/regions/{{figure}}/libsizenorm/{condition}-v-{control}/{{figure}}_{condition}-v-{control}_libsizenorm-individual-occupancy-heatmaps.svg", zip, condition=conditiongroups, control=controlgroups), figure=QUANT) if sipassing else expand(expand("nucleosome_calling/regions/{{figure}}/libsizenorm/{condition}-v-{control}/{{figure}}_{condition}-v-{control}_libsizenorm-individual-occupancy-heatmaps.svg", zip, condition=conditiongroups, control=controlgroups), figure=QUANT),
        ##differential nucleosome levels over transcripts
        #expand("diff_levels/{condition}-v-{control}/libsizenorm/{condition}-v-{control}-results-libsizenorm-all.tsv", zip, condition=conditiongroups, control=controlgroups) + expand("diff_levels/{condition}-v-{control}/spikenorm/{condition}-v-{control}-results-spikenorm-all.tsv", zip, condition=conditiongroups_si, control=controlgroups_si) if sipassing else expand("diff_levels/{condition}-v-{control}/libsizenorm/{condition}-v-{control}-results-libsizenorm-all.tsv", zip, condition=conditiongroups, control=controlgroups)

def plotcorrsamples(wc):
    if wc.condition=="all":
        if wc.norm=="libsizenorm": #condition==all,norm==lib
            return list(SAMPLES.keys())
        else: #condition==all,norm==spike
            return list(sisamples.keys())
    elif wc.norm=="libsizenorm": #condition!=all;norm==lib
        return [k for k,v in PASSING.items() if v["group"] in (wc.control, wc.condition)]
    else: #condition!=all;norm==spike
        return [k for k,v in sipassing.items() if v["group"] in (wc.control, wc.condition)]

def cluster_samples(status, norm, cluster_groups):
    dd = SAMPLES if status=="all" else PASSING
    if norm=="libsizenorm": #condition!=all;norm==lib
        return [k for k,v in dd.items() if v["group"] in cluster_groups]
    else: #condition!=all;norm==spike
        return [k for k,v in dd.items() if v["group"] in cluster_groups and v["spikein"]=="y"]

def getsamples(ctrl, cond):
    return [k for k,v in PASSING.items() if (v["group"]==ctrl or v["group"]==cond)]

include: "rules/mnase-seq_clean_reads.smk"
include: "rules/mnase-seq_alignment.smk"
include: "rules/mnase-seq_fastqc.smk"
include: "rules/mnase-seq_library_processing_summary.smk"
include: "rules/mnase-seq_genome_coverage.smk"
include: "rules/mnase-seq_sample_similarity.smk"
include: "rules/mnase-seq_datavis.smk"
include: "rules/mnase-seq_quantification.smk"
include: "rules/mnase-seq_differential_occupancy.smk"

