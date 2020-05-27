#!/usr/env/python

import argparse
import numpy as np
import pybedtools as pybt
import pyBigWig as pybw

def transcript_nucleosome_intersect(transcript_annotation_path,
                                    nucleosome_annotation_path):
    return (pybt.BedTool(transcript_annotation_path)
            .intersect(pybt.BedTool(nucleosome_annotation_path).remove_invalid(),
                       wo=True,
                       F=0.5)
            .to_dataframe(names=["chrom",
                                 "transcript_start",
                                 "transcript_end",
                                 "transcript_name",
                                 "score",
                                 "strand",
                                 "nuc_chrom",
                                 "nuc_start",
                                 "nuc_end",
                                 "nuc_dyad_position",
                                 "nuc_dyad_value",
                                 "nuc_fuzziness",
                                 "overlap"])[["chrom",
                                              "transcript_start",
                                              "transcript_end",
                                              "transcript_name",
                                              "strand",
                                              "nuc_start",
                                              "nuc_end"]])

def get_nucleosome_mode(coverage_objects,
                        nucleosome_feature):
    normalized_coverages = []
    for coverage_object in coverage_objects:
        normalized_coverages.append(coverage_object.values(nucleosome_feature["chrom"],
                                                           nucleosome_feature["nuc_start"],
                                                           nucleosome_feature["nuc_end"],
                                                           numpy=True))
    mean_normalized_coverage = np.mean(normalized_coverages, axis=0)
    return nucleosome_feature["nuc_start"] + np.argmax(mean_normalized_coverage)

def augment_distance_information(df,
                                coverage_objects,
                                output_tsv_path):
    df["nuc_mode"] = df.apply(lambda x: get_nucleosome_mode(coverage_objects, x),
                              axis=1)
    df["distance_from_TSS"] = df.apply(lambda x: x["nuc_mode"] - x["transcript_start"] if x["strand"]=="+"
                                       else x["transcript_end"] - (x["nuc_mode"] + 1),
                                       axis=1)
    df = df[df["distance_from_TSS"] >= 0]
    df = df.sort_values(by=["chrom",
                            "transcript_start",
                            "distance_from_TSS"])
    df["nuc_index"] = df.groupby(["chrom",
                                  "transcript_start",
                                  "transcript_end",
                                  "transcript_name",
                                  "strand"]).cumcount()
    df["distance_from_previous"] = (df["distance_from_TSS"] -
                                    (df.groupby(["chrom",
                                                "transcript_start",
                                                "transcript_end",
                                                "transcript_name",
                                                "strand"])
                                     .shift(1,
                                            fill_value=0))["distance_from_TSS"])
    df.to_csv(output_tsv_path,
              index=False,
              sep="\t")
    return df

def main(transcript_annotation_path = "Scer_polIItranscripts-adjustedTSS.bed",
         nucleosome_annotation_path =
         "nucleosome_quantification_data_WT-37C.Fnor.smooth.positions.xls",
         coverage_paths =
         ["WT-37C-1_mnase-midpoint_smoothed-libsizenorm.bw",
          "WT-37C-2_mnase-midpoint_smoothed-libsizenorm.bw"],
        output_tsv_path="control.tsv"):

    df = transcript_nucleosome_intersect(
        transcript_annotation_path=transcript_annotation_path,
        nucleosome_annotation_path=nucleosome_annotation_path)

    coverage_objects = [pybw.open(coverage_path) for coverage_path in
            coverage_paths]

    df = augment_distance_information(df=df,
            coverage_objects=coverage_objects,
            output_tsv_path=output_tsv_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get nucleosome distances relative to TSSs, given transcript annotation, DANPOS nucleosome calls, and smoothed MNase dyad coverage.")
    parser.add_argument('-t', dest='transcript_annotation', type=str, \
            help='Path to transcript annotation, in BED6 format.')
    parser.add_argument('-n', dest='nucleosome_annotation', type=str, \
            help='Path to DANPOS nucleosome calls (TSV labeled as XLS).')
    parser.add_argument('-c', dest='coverage', type=str, nargs="+", \
            help='Paths to smoothed MNase dyad coverage files, in BigWig format.')
    parser.add_argument('-o', dest='output_tsv', type=str, \
            help='Path to output TSV file.')
    args = parser.parse_args()
    main(transcript_annotation_path = args.transcript_annotation ,
         nucleosome_annotation_path = args.nucleosome_annotation,
         coverage_paths = args.coverage,
         output_tsv_path = args.output_tsv)

