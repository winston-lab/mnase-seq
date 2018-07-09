#!/usr/bin/env python

import sys
import csv

# if using spikein normalization, set control to 10M reads and condition to 10M*(condition SI pct)/(control SI pct)
def danpos_norm(norm, condition, control,  condition_samples, control_samples, si_table):
    if norm=="spikenorm":
        cond_count, ctrl_count, cond_val, ctrl_val = 0,0,0,0
        with open(si_table) as si_table:
            si_table = csv.reader(si_table, delimiter="\t")
            for row in si_table:
                if row[0] in condition_samples.split(","):
                    cond_val += float(row[4])/float(row[2])
                    cond_count += 1
                if row[0] in control_samples.split(","):
                    ctrl_val += float(row[4])/float(row[2])
                    ctrl_count += 1
        cond_sipct = cond_val/cond_count
        ctrl_sipct = ctrl_val/ctrl_count
        spikein_counts = int(1e7*ctrl_sipct*(1-cond_sipct)/((1-ctrl_sipct)*cond_sipct))
        spikein_string = "--count nucleosome_quantification/data/{condition}/:{spikein_counts},nucleosome_quantification/data/{control}/:1e7".format(condition=condition, spikein_counts=spikein_counts, control=control)
        return spikein_string
    else:
        return ""

print(danpos_norm(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]))

