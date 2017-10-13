#!/usr/bin/env python

import argparse
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gsmooth
import pyBigWig as pybw

parser = argparse.ArgumentParser(description='Smooth bigwig file with Gaussian kernel of given bandwidth.')
parser.add_argument('-b', dest = 'bandwidth', type=int, default = 20, help='Gaussian kernel bandwidth (standard deviation)')
parser.add_argument('-i', dest = 'infile', type=str, help='path to input BigWig')
parser.add_argument('-o', dest = 'outfile', type=str, help='path to smoothed output BigWig')
args = parser.parse_args()

inbw = pybw.open(args.infile)
outbw = pybw.open(args.outfile, "w")

outbw.addHeader(list(inbw.chroms().items()))

for chrom in inbw.chroms():
    raw = inbw.values(chrom, 0, inbw.chroms(chrom), numpy=True)
    smoothed = gsmooth(raw, sigma=args.bandwidth, order=0, mode='mirror')
    outbw.addEntries(chrom, 0, values=smoothed, span=1, step=1)

inbw.close()
outbw.close()
