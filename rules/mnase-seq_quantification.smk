#!/usr/bin/env python

localrules:
    group_bam_for_danpos,
    danpos_over_annotations,
    cat_danpos_over_annotations

#danpos has an unfriendly syntax for inputs, requires putting replicates into same directory
rule group_bam_for_danpos:
    input:
        "alignment/{sample}_mnase-seq-experimental.bam" if SISAMPLES else "alignment/{sample}_mnase-seq.bam",
    output:
        temp("nucleosome_quantification/data/{group}/{sample}.bam")
    shell: """
        cp {input} {output}
        """

# if using spikein normalization, set control to 10M reads and condition to 10M*(condition SI pct)/(control SI pct)
def danpos_norm(norm, condition, control, si_table):
    if norm=="spikenorm":
        cond_count, ctrl_count, cond_val, ctrl_val = 0,0,0,0
        with open(si_table) as si_table:
            si_table = csv.reader(si_table, delimiter="\t")
            for row in si_table:
                if row[0] in [k for k,v in SIPASSING.items() if v["group"]==condition]:
                    # vals = [int(x) for x in row[2].split()]
                    #TODO: should really fix the sicounts file to be a proper tsv file...
                    # cond_val += vals[2]/vals[0]
                    cond_val += float(row[4])/float(row[2])
                    cond_count += 1
                if row[0] in [k for k,v in SIPASSING.items() if v["group"]==control]:
                    # vals = [int(x) for x in row[2].split()]
                    # ctrl_val += vals[2]/vals[0]
                    ctrl_val += float(row[4])/float(row[2])
                    ctrl_count += 1
        cond_sipct = cond_val/cond_count
        ctrl_sipct = ctrl_val/ctrl_count
        spikein_counts = int(1e7*ctrl_sipct*(1-cond_sipct)/((1-ctrl_sipct)*cond_sipct))
        spikein_string = f"--count nucleosome_quantification/data/{condition}/:{spikein_counts},nucleosome_quantification/data/{control}/:1e7"
        return spikein_string
    else:
        return ""

rule danpos_quantification:
    input:
        bam = lambda wc: ["nucleosome_quantification/data/" + PASSING[x]['group'] + "/" + x + ".bam" for x in PASSING] if wc.norm=="libsizenorm" else ["nucleosome_quantification/data/" + SIPASSING[x]['group'] + "/" + x + ".bam" for x in SIPASSING],
        si_table = lambda wc: [] if wc.norm=="libsizenorm" else "qual_ctrl/spikein/mnase-seq_spikein-counts.tsv"
    output:
        "nucleosome_quantification/{condition}-v-{control}/{norm}/nucleosome_quantification_data_{condition}-nucleosome_quantification_data_{control}.positions.integrative.xls",
        "nucleosome_quantification/{condition}-v-{control}/{norm}/reference_positions.xls",
        "nucleosome_quantification/{condition}-v-{control}/{norm}/diff/nucleosome_quantification_data_{condition}-nucleosome_quantification_data_{control}.pois_diff.wig",
        "nucleosome_quantification/{condition}-v-{control}/{norm}/pooled/nucleosome_quantification_data_{condition}.Fnor.smooth.positions.ref_adjust.xls",
        "nucleosome_quantification/{condition}-v-{control}/{norm}/pooled/nucleosome_quantification_data_{condition}.Fnor.smooth.positions.xls",
        "nucleosome_quantification/{condition}-v-{control}/{norm}/pooled/nucleosome_quantification_data_{condition}.Fnor.smooth.wig",
        "nucleosome_quantification/{condition}-v-{control}/{norm}/pooled/nucleosome_quantification_data_{control}.Fnor.smooth.positions.ref_adjust.xls",
        "nucleosome_quantification/{condition}-v-{control}/{norm}/pooled/nucleosome_quantification_data_{control}.Fnor.smooth.positions.xls",
        "nucleosome_quantification/{condition}-v-{control}/{norm}/pooled/nucleosome_quantification_data_{control}.Fnor.smooth.wig"
    params:
        spikein_string = lambda wc: danpos_norm(wc.norm, wc.condition, wc.control, "qual_ctrl/spikein/mnase-seq_spikein-counts.tsv")
    conda:
        "../envs/danpos.yaml"
    log: "logs/danpos_quantification/danpos_quantification-{condition}-v-{control}-{norm}.log"
    shell: """
        (python scripts/danpos-2.2.2/danpos.py dpos nucleosome_quantification/data/{wildcards.condition}/:nucleosome_quantification/data/{wildcards.control}/ {params.spikein_string} --paired 1 --edge 1 --span 1 -o nucleosome_quantification/{wildcards.condition}-v-{wildcards.control}/{wildcards.norm}) &> {log}
        """

rule plot_danpos_results:
    input:
        results = "nucleosome_quantification/{condition}-v-{control}/{norm}/nucleosome_quantification_data_{condition}-nucleosome_quantification_data_{control}.positions.integrative.xls",
    output:
        shift_hist = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-dyad-shift-histogram.svg",
        occupancy_volcano = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-occupancy_volcano.svg",
        fuzziness_volcano = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-fuzziness_volcano.svg",
        occupancy_v_fuzziness = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-occupancy-v-fuzziness-scatter.svg",
        occupancy_violin = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-occupancy-violin.svg",
        occupancy_freqpoly = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-occupancy-freqpoly.svg",
        occupancy_ecdf = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-occupancy-ecdf.svg",
        fuzziness_violin = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-fuzziness-violin.svg",
        fuzziness_freqpoly = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-fuzziness-freqpoly.svg",
        fuzziness_ecdf = "nucleosome_quantification/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-fuzziness-ecdf.svg"
    conda: "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_danpos_results.R"

rule danpos_over_annotations:
    input:
        indiv_condition = "nucleosome_quantification/{condition}-v-{control}/{norm}/pooled/nucleosome_quantification_data_{condition}.Fnor.smooth.positions.xls",
        indiv_control = "nucleosome_quantification/{condition}-v-{control}/{norm}/pooled/nucleosome_quantification_data_{control}.Fnor.smooth.positions.xls",
        integrated = "nucleosome_quantification/{condition}-v-{control}/{norm}/nucleosome_quantification_data_{condition}-nucleosome_quantification_data_{control}.positions.integrative.xls",
        annotation = lambda wc: QUANT[wc.figure]["annotations"][wc.annotation]["path"],
        fasta = config["genome"]["fasta"]
    output:
        individual = temp("nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{annotation}_{norm}-individual.tsv"),
        integrated = temp("nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{annotation}_{norm}-integrated.tsv")
    params:
        upstream = lambda wc: QUANT[wc.figure]["upstream"],
        dnstream = lambda wc: QUANT[wc.figure]["dnstream"],
        label = lambda wc: QUANT[wc.figure]["annotations"][wc.annotation]["label"]
    log: "logs/danpos_over_annotations/danpos_over_annotations_{condition}-v-{control}_{figure}-{annotation}-{norm}.log"
    shell: """
        (bedtools slop -i {input.annotation} -g <(faidx {input.fasta} -i chromsizes) -l {params.upstream} -r {params.dnstream} -s | bedtools intersect -a stdin -b <(awk 'BEGIN{{FS=OFS="\t"}} NR>1 && $2>0' {input.indiv_control}) -wo | awk 'BEGIN{{FS=OFS="\t"}}{{print $0, "{wildcards.control}", "{params.label}"}}' | \
        cat - <(bedtools slop -i {input.annotation} -g <(faidx {input.fasta} -i chromsizes) -l {params.upstream} -r {params.dnstream} -s | bedtools intersect -a stdin -b <(awk 'BEGIN{{FS=OFS="\t"}} NR>1 && $2>0' {input.indiv_condition}) -wo | awk 'BEGIN{{FS=OFS="\t"}}{{print $0, "{wildcards.condition}", "{params.label}"}}') > {output.individual}) &> {log}
        (bedtools slop -i {input.annotation} -g <(faidx {input.fasta} -i chromsizes) -l {params.upstream} -r {params.dnstream} -s | bedtools intersect -a stdin -b <(awk 'BEGIN{{FS=OFS="\t"}} NR>1 && $2>0' {input.integrated}) -wo | awk 'BEGIN{{FS=OFS="\t"}}{{print $0, "{params.label}"}}'> {output.integrated}) &>> {log}
        """

rule cat_danpos_over_annotations:
    input:
        individual = lambda wc: expand("nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{annotation}_{norm}-individual.tsv", annotation=[k for k,v in QUANT[wc.figure]["annotations"].items()], condition=wc.condition, control=wc.control, figure=wc.figure, norm=wc.norm),
        integrated = lambda wc: expand("nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{annotation}_{norm}-integrated.tsv", annotation=[k for k,v in QUANT[wc.figure]["annotations"].items()], condition=wc.condition, control=wc.control, figure=wc.figure, norm=wc.norm)
    output:
        individual = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}-allannotations-individual-{norm}.tsv",
        integrated = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}-allannotations-integrated-{norm}.tsv",
    log: "logs/cat_danpos_over_annotations/cat_danpos_over_annotations_{condition}-v-{control}_{figure}-{norm}.log"
    shell: """
        (cat <(echo -e "chrom\tfeat_start\tfeat_end\tfeat_name\tfeat_score\tfeat_strand\tnuc_start\tnuc_end\tnuc_summit\toccupancy\tfuzziness\tgroup\tannotation") <(cut -f7,13 --complement {input.individual}) > {output.individual}) &> {log}
        (cat <(echo -e "feat_chrom\tfeat_start\tfeat_end\tfeat_name\tfeat_score\tfeat_strand\tnuc_chrom\tnuc_start\tnuc_end\tnuc_center\tctrl_summit_loc\tcond_summit_loc\tdiff_summit_loc\tcond_ctrl_dist\tctrl_summit_val\tcond_summit_val\tsummit_lfc\tsummit_diff_logpval\tsummit_diff_fdr\tctrl_point_val\tcond_point_val\tpoint_lfc\tpoint_diff_logpval\tpoint_diff_fdr\tctrl_fuzziness\tcond_fuzziness\tfuzziness_lfc\tfuzziness_diff_logpval\tfuzziness_diff_fdr\toverlap\tannotation") {input.integrated} >  {output.integrated}) &> {log}
        """

rule danpos_vis_over_annotations:
    input:
        individual = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}-allannotations-individual-{norm}.tsv",
        integrated = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}-allannotations-integrated-{norm}.tsv",
        annotations = lambda wc: [v["path"] for k,v in QUANT[wc.figure]["annotations"].items()]
    output:
        indiv_occ_hmap = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-individual-occupancy-heatmaps.svg",
        indiv_fuzz_hmap = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-individual-fuzziness-heatmaps.svg",
        indiv_occ_meta = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-individual-occupancy-metagene.svg",
        indiv_fuzz_meta = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-individual-fuzziness-metagene.svg",
        integrated_occ_summit_hmap = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-integrated-summit-occupancy-heatmap.svg",
        integrated_occ_point_hmap = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-integrated-point-occupancy-heatmap.svg",
        integrated_fuzz_hmap = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-integrated-fuzziness-heatmap.svg",
        integrated_displacement_hmap = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-integrated-displacement-heatmap.svg",
        integrated_displacement_segment_hmap = "nucleosome_quantification/regions/{norm}/{figure}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-integrated-displacement-segment-heatmap.svg",
        integrated_occ_summit_meta = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-integrated-summit-occupancy-metagene.svg",
        integrated_occ_point_meta = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-integrated-point-occupancy-metagene.svg",
        integrated_fuzz_meta = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-integrated-fuzziness-metagene.svg",
        integrated_displacement_meta = "nucleosome_quantification/regions/{figure}/{norm}/{condition}-v-{control}/{figure}_{condition}-v-{control}_{norm}-integrated-displacement-metagene.svg",
    params:
        anno_labels = lambda wc: [v["label"] for k,v in QUANT[wc.figure]["annotations"].items()],
        refpoint = lambda wc: QUANT[wc.figure]["refpoint"],
        refptlabel = lambda wc: QUANT[wc.figure]["refpointlabel"],
        sortmethod = lambda wc: QUANT[wc.figure]["arrange"],
        binsize = lambda wc: QUANT[wc.figure]["binsize"],
        occupancy_cutoffs = lambda wc: [QUANT[wc.figure]["occupancy_cutoff_low"], QUANT[wc.figure]["occupancy_cutoff_high"]],
        fuzziness_cutoffs = lambda wc: [QUANT[wc.figure]["fuzziness_cutoff_low"], QUANT[wc.figure]["fuzziness_cutoff_high"]],
        occupancy_lfc_limit = lambda wc: QUANT[wc.figure]["occupancy_lfc_limit"],
        fuzziness_lfc_limit = lambda wc: QUANT[wc.figure]["fuzziness_lfc_limit"],
        displacement_limit = lambda wc: QUANT[wc.figure]["displacement_limit"],
        upstream = lambda wc: QUANT[wc.figure]["upstream"],
        max_length = lambda wc: QUANT[wc.figure]["max_length"],
        trim_pct = lambda wc: QUANT[wc.figure]["trim_pct"],
    conda: "../envs/tidyverse.yaml"
    script:
        "../scripts/mnase_quant_vis.R"

