# -*- coding: utf-8 -*-
"""
@author: Luis Quintales - IBFG - Universidad de Salamanca
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

def parsing () :
  import argparse
  import os
  parser = argparse.ArgumentParser("python nucwave_pe.py")
  parser.add_argument("-w", "--wigfiles", help="write intermediate wig files",action="store_true")
  parser.add_argument("-o", metavar="OUTDIR", help="output directory",default=".")
  parser.add_argument("-g", metavar="GENOMEFILE", help="FASTA genome file",
                      required=True)
  parser.add_argument("-a", metavar="ALIGNFILE", help="BOWTIE alignment file",
                      required=True)   
  parser.add_argument("-p", metavar="PREFIXNAME", help="Prefix name for result files",
                      required=True)                     
  args=parser.parse_args()
  if not os.path.isdir(args.o) :  
    try :
      os.mkdir(args.o) 
    except OSError as e:
     print('Making output directory "%s": %s' % ( args.o , e[1]))
     exit()
  if not os.access(args.o,os.W_OK) :
    print('I can not write in output directory "%s"' % (args.o))
    exit()
  if not os.path.isfile(args.g) :
    print('Genomic file "%s" does not exist' % (args.g))
    exit()  
  if not os.access(args.g,os.R_OK) :
    print('No permission to read genomic file "%s"' % (args.g))
    exit()
  if not os.path.isfile(args.a) :
    print('Alignment file "%s" does not exist' % (args.a))
    exit()  
  if not os.access(args.a,os.R_OK) :
    print('No permission to read aligment file "%s"' % (args.a))
    exit()    

  return args.o, args.g, args.a, args.p, args.wigfiles

def buscapicos (signals,umbral,vecindad) :
  picos={} 
  cont_picos = 0 
  pend_ant = 1 if signals[0] < signals[1] else -1
  for coord in range(1,len(signals)-1) :
    pend = 1 if signals[coord] < signals[coord+1] else -1
    if pend_ant == 1 and pend == -1 and signals[coord] > umbral:
      picos[coord] = signals[coord]
      cont_picos += 1 
    pend_ant = pend
  
  picos_ok = {} 
  picos_ord =  sorted(picos.keys())
  for i in range(1,len(picos_ord)-1) :
    if (picos_ord[i]-picos_ord[i-1] > vecindad and \
        picos_ord[i+1]-picos_ord[i] > vecindad) :
      coord = picos_ord[i]
      picos_ok[coord] = signals[coord]
  return picos_ok

def denoising (signal) :
  level = 5
  coeffs = pywt.wavedec(signal, 'bior3.1',level=level)
  for i in range(1,level+1) : 
    #coeffs[i] = [0]*len(coeffs[i])
    coeffs[i] = np.zeros_like(coeffs[i])
  signal_wave = pywt.waverec(coeffs, 'bior3.1')
  ls = len(signal)
  lsw = len(signal_wave)
  if ls < lsw :  
    signal_wave = np.delete(signal_wave,-1)
  return signal_wave
  
def fasta_chr (fasta_file) :
  chrseq = {} 
  chrname = "" 
  fichero = open(fasta_file)
  primero = True 
  for linea in fichero :
    linea = linea.strip("\n") 
    if linea[0] == ">" :
      if primero :
        primero = False
      else :
        chrseq[chrname] =  sequence.replace(" ", "").upper()
      chrname = linea.replace(">","").split()[0]
      sequence=""  
    else :
      sequence += linea
  chrseq[chrname] =  sequence.replace(" ", "").upper()
  return chrseq
  
### MAIN PROGRAM
import time  
import numpy as np
import pywt

# Parsing command line
(OUTdir,FASTA_file,BOWTIE_file,exper,wig_int) = parsing()
WIG_file = '%s/%s' % (OUTdir, exper)

minsize = 40
maxsize = 200
extension = minsize//2

# 
time0 = int(time.time())
print("Reading and processing genome")
chrseq = fasta_chr (FASTA_file)
chrlen = {}
for chrid in chrseq.keys() :
  chrlen[chrid] = len(chrseq[chrid])
print("\tTime taken: %d seconds" % (time.time()-time0))

# 
time0 = int(time.time())
print("Reading and processing alignment file")
sizesreads = [0]*maxsize
iniciosPE = {} 
finalesPE = {} 
centrosPE = {} 
iniciosPEext = {} 
finalesPEext = {} 
depth_trimming_PE = {} 
for chrid in chrseq.keys() :
  iniciosPE[chrid] = [0]*chrlen[chrid]
  finalesPE[chrid] = [0]*chrlen[chrid]
  centrosPE[chrid] = [0]*chrlen[chrid]
  iniciosPEext[chrid] = [0]*chrlen[chrid]
  finalesPEext[chrid] = [0]*chrlen[chrid]
  depth_trimming_PE[chrid] = [0]*chrlen[chrid]
  
# 
F = open (BOWTIE_file)
nfragmentos=0
nfragmentos_filtrados=0
for linea1 in F :
  strand1, chrid1, coord1, seq1 = linea1.strip('\n').split('\t')
  linea2 = next(F)
  strand2, chrid2, coord2, seq2 = linea2.strip('\n').split('\t')
  coord1 = int(coord1)
  coord2 = int(coord2)
  nfragmentos += 1

  iniPE = coord1
  finPE = coord2 + len(seq2) - 1
  cenPE = (iniPE+finPE)//2

  d2d_long = finPE - iniPE + 1
  if d2d_long < maxsize and d2d_long > minsize:     
#  if d2d_long < maxsize :     
    nfragmentos_filtrados += 1
    sizesreads[d2d_long] += 1
    iniciosPE[chrid1][iniPE] += 1
    finalesPE[chrid1][finPE] += 1
    centrosPE[chrid1][cenPE] += 1
    iniciosPEext[chrid1][cenPE-extension] += 1
    finalesPEext[chrid1][cenPE+extension] += 1
F.close()
print("\tNumber of reads processed: %d" % nfragmentos_filtrados)
print("\tTime taken: %d seconds" % (time.time()-time0))
time0 = int(time.time())

#
if wig_int:
  time0 = int(time.time())
  print("Writing intermediate files: fragment size histogram")
  suffix = 'historeadsize'
  SIZES = open(WIG_file + ".%s" % suffix,"w")
  for i in range(len(sizesreads)) :
     SIZES.write('%d\t%d\n' % (i+1,sizesreads[i]) ) 
  SIZES.close()
  print("\tTime taken: %d seconds" % (time.time()-time0))

#
if wig_int:
  time0 = int(time.time())
  print("Writing intermediate files: cut points per strand")

  suffix = 'cut_p'
  WIG = open(WIG_file + "_%s.wig" % suffix,"w")
  for chrid in chrseq.keys() :
    WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    for i in range(chrlen[chrid]) :
      WIG.write('%d\n' % iniciosPE[chrid][i] )
  WIG.close()
  
  suffix = 'cut_m'
  WIG = open(WIG_file + "_%s.wig" % suffix,"w")
  for chrid in chrseq.keys() :
    WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    for i in range(chrlen[chrid]) :
      WIG.write('%d\n' % finalesPE[chrid][i] )
  WIG.close()
  print("\tTime taken: %d seconds" % (time.time()-time0))

#
if wig_int:
  time0 = int(time.time())
  print("Writing intermediate files: depth coverage for complete PE reads")
  suffix = 'depth_complete_PE'
  WIG = open(WIG_file + "_%s.wig" % suffix,"w")
  for chrid in chrseq.keys() :
    WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    deep = 0
    for i in range(chrlen[chrid]) :
      deep = deep+iniciosPE[chrid][i]-finalesPE[chrid][i]
      WIG.write('%d\n' % deep )
  WIG.close()
  print("\tTime taken: %d seconds" % (time.time()-time0))

#
if wig_int:
  time0 = int(time.time())
  print("Writing intermediate files: PE center count")
  suffix = 'PEcenter'
  WIG = open(WIG_file + "_%s.wig" % suffix,"w")
  for chrid in chrseq.keys() :
    WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    for i in range(chrlen[chrid]) :
      WIG.write('%d\n' % centrosPE[chrid][i] )
  WIG.close()
  print("\tTime taken: %d seconds" % (time.time()-time0))

#
time0 = int(time.time())
print("Depth coverage for trimmed PE reads")
if wig_int:
  print("Writing intermediate files: depth coverage for trimmed PE reads")

suffix = 'depth_trimmed_PE'
if wig_int: WIG = open(WIG_file + "_%s.wig" % suffix,"w")
for chrid in chrseq.keys() :
  if wig_int: WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
  if wig_int: WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
  deep = 0
  for i in range(chrlen[chrid]) :
    deep = deep+iniciosPEext[chrid][i]-finalesPEext[chrid][i]
    depth_trimming_PE[chrid][i] = deep
    if wig_int: WIG.write('%d\n' % deep )
if wig_int: WIG.close()
print("\tTime taken: %d seconds" % (time.time()-time0))

# Wavelet denoising 
time0 = int(time.time())
print("Depth coverage wavelet denoising")
depth_wl_trimming_PE = {}
for chrid in chrseq.keys() :
  depth_wl_trimming_PE[chrid] = denoising(depth_trimming_PE[chrid])
print("\tTime taken: %d seconds" % (time.time()-time0))

time0 = int(time.time())
print("Mean depth coverage calculation")
suma = 0
count = 0
for chrid in chrseq.keys() :
  for i in range(len(depth_wl_trimming_PE[chrid])) :
    if depth_wl_trimming_PE[chrid][i] <= 0.0 : 
      depth_wl_trimming_PE[chrid][i] = 0.0 
    else :  
      count += 1
      suma += depth_wl_trimming_PE[chrid][i] 
mean_depth = suma/count 
print("\tMean depth coverage: %.0f" % mean_depth)
print("\tTime taken: %d seconds" % (time.time()-time0))     
# 
time0 = int(time.time())
print("Writing results file: Depth coverage wavelet denoised and normalized")
suffix = 'depth_wl_trimmed_PE'
WIG = open(WIG_file + "_%s.wig" % suffix,"w")
for chrid in chrseq.keys() :
  WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
  WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
  for i in range(chrlen[chrid]) :
    deep_n = depth_wl_trimming_PE[chrid][i]/mean_depth
    WIG.write('%.4f\n' % deep_n )
WIG.close()
print("\tTime taken: %d seconds" % (time.time()-time0))








