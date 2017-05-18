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
  parser = argparse.ArgumentParser("python nucwave_sr.py")
  parser.add_argument("-w", "--wigfiles", help="write intermediate wig files",action="store_true")
  parser.add_argument("-c", "--chemical", help="chemical cleavage",action="store_true")
  parser.add_argument("-o", metavar="OUTDIR", help="output directory",default=".")
  parser.add_argument("-g", metavar="GENOMEFILE", help="FASTA genome file",
                      required=True)
  parser.add_argument("-a", metavar="ALIGNFILE", help="BOWTIE alignment file",
                      required=True)   
  parser.add_argument("-p", metavar="PREFIXNAME", help="Prefix name for result files",
                      required=True)                     
  args=parser.parse_args()
  ## OUTPUT DIRECTORY
  # Create if not exists
  if not os.path.isdir(args.o) :  
    try :
      os.mkdir(args.o) 
    except OSError as e:
     print('Making output directory "%s": %s' % ( args.o , e[1]))
     exit()
  # Verify write permission in output diectory
  if not os.access(args.o,os.W_OK) :
    print('I can not write in output directory "%s"' % (args.o))
    exit()
  # Verify genome file existence and access
  if not os.path.isfile(args.g) :
    print('Genomic file "%s" does not exist' % (args.g))
    exit()  
  if not os.access(args.g,os.R_OK) :
    print('No permission to read genomic file "%s"' % (args.g))
    exit()
  # Verify genome file existence and access
  if not os.path.isfile(args.a) :
    print('Alignment file "%s" does not exist' % (args.a))
    exit()  
  if not os.access(args.a,os.R_OK) :
    print('No permission to read aligment file "%s"' % (args.a))
    exit()    

  return args.o, args.g, args.a, args.p, args.wigfiles, args.chemical

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
(OUTdir,FASTA_file,BOWTIE_file,exper,wig_int,chemical) = parsing()
WIG_file = '%s/%s' % (OUTdir, exper)

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
fir = {} # First position in read
lir = {} # Last position in read
for chrid in chrseq.keys() :
  fir[chrid] = {}
  fir[chrid]['+'] = [0]*chrlen[chrid]
  fir[chrid]['-'] = [0]*chrlen[chrid]
  lir[chrid] = {}
  lir[chrid]['+'] = [0]*chrlen[chrid]
  lir[chrid]['-'] = [0]*chrlen[chrid]

# 
F = open (BOWTIE_file)
nreads=0
switch_strand = {'+': '-','-' :'+'}
for linea in F :
  strand, chrid, coord, seq = linea.strip('\n').split('\t')
  if chemical :
    strand = switch_strand[strand]
  nreads += 1
  coord = int(coord)
  fir[chrid][strand][coord] += 1
  lir[chrid][strand][coord+len(seq)-1] += 1
F.close()
print("\tNumber of reads processed: %d" % nreads)
print( "\tTime taken: %d seconds" % (time.time()-time0))

if wig_int:
  time0 = int(time.time())
  print("Writing intermediate files: cut points per strand")
  suffix = 'cut_p'
  WIG = open(WIG_file + "_%s.wig" % suffix,"w")
  for chrid in chrseq.keys() :
    WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    for i in range(chrlen[chrid]) :
      WIG.write('%d\n' % fir[chrid]['+'][i] )
  WIG.close()

  suffix = 'cut_m'
  WIG = open(WIG_file + "_%s.wig" % suffix,"w")
  for chrid in chrseq.keys() :
    WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    for i in range(chrlen[chrid]) :
      WIG.write('%d\n' % lir[chrid]['-'][i])
  WIG.close()
  print("\tTime taken: %d seconds" % (time.time()-time0))

# 
time0 = int(time.time())
print("Depth coverage calculation")
depth = {} # 
depth_wl = {} # 
for chrid in chrseq.keys() :
  # 
  depth[chrid] = {} 
  depth[chrid]['+'] = []
  depth[chrid]['-'] = []
  depth_wl[chrid] = {} 
  depth_wl[chrid]['+'] = []
  depth_wl[chrid]['-'] = []
  # 
  depth_pos = 0 # 
  for i in range(chrlen[chrid]) :
    depth_pos = depth_pos+fir[chrid]['+'][i]-lir[chrid]['+'][i]
    depth[chrid]['+'].append(depth_pos)
  depth_pos = 0 # 
  for i in range(chrlen[chrid]) :
    depth_pos = depth_pos+fir[chrid]['-'][i]-lir[chrid]['-'][i]
    depth[chrid]['-'].append(depth_pos)
print("\tTime taken: %d seconds" % (time.time()-time0))    

# 
time0 = int(time.time())
print("Depth coverage denoising")
for chrid in chrseq.keys() : 
  depth_wl[chrid]['+'] = denoising(depth[chrid]['+'])
  depth_wl[chrid]['-'] = denoising(depth[chrid]['-'])
print("\tTime taken: %d seconds" % (time.time()-time0))

if wig_int:
  time0 = int(time.time())
  print("Writing intermediate files: depth coverage per strand")
  suffix = 'depth_p'
  WIG = open(WIG_file + "_%s.wig" % suffix,"w")
  for chrid in chrseq.keys() :
    WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    for d in depth[chrid]['+'] :
      WIG.write('%d\n' % d )
  WIG.close()
  
  suffix = 'depth_m'
  WIG = open(WIG_file + "_%s.wig" % suffix,"w")
  for chrid in chrseq.keys() :
    WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    for d in depth[chrid]['-'] :
      WIG.write('%d\n' % d )
  WIG.close()
  print("\tTime taken: %d seconds" % (time.time()-time0))

time0 = int(time.time())
print("Mean depth coverage calculation")
if wig_int:
  print("Writing intermediate files: Depth coverage wavelet denoised per strand")

suma = 0  
count = 0 
suffix = 'depth_p_wl'
if wig_int: WIG = open(WIG_file + "_%s.wig" % suffix,"w")
for chrid in chrseq.keys() :
  if wig_int: WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
  if wig_int: WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
  for i in range(len(depth_wl[chrid]['+'])) :
    if depth_wl[chrid]['+'][i] <= 0 : 
      depth_wl[chrid]['+'][i] = 0
    else : 
      count += 1
      suma += depth_wl[chrid]['+'][i]
    if wig_int: WIG.write('%.4f\n' % depth_wl[chrid]['+'][i] )
if wig_int: WIG.close()

suffix = 'depth_m_wl'
if wig_int: WIG = open(WIG_file + "_%s.wig" % suffix,"w")
for chrid in chrseq.keys() :
  if wig_int: WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
  if wig_int: WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
  for i in range(len(depth_wl[chrid]['-'])) :
    if depth_wl[chrid]['-'][i] < 0 : 
      depth_wl[chrid]['-'][i] = 0
    else : 
      count += 1
      suma += depth_wl[chrid]['-'][i]
    if wig_int: WIG.write('%.4f\n' % depth_wl[chrid]['-'][i] )
if wig_int: WIG.close()

mean_wl = suma/count
print("\tMean depth coverage: %.0f" % mean_wl)
print("\tTime taken: %d seconds" % (time.time()-time0))

#
time0 = int(time.time())
print("Shift calculation")
suma = 0
count = 0
for chrid in chrseq.keys() :
  picos = {}
  picos['+'] = buscapicos(depth_wl[chrid]['+'], 2.0*mean_wl, 100)
  picos['-'] = buscapicos(depth_wl[chrid]['-'], 2.0*mean_wl, 100)
  for coord in picos['+'] :
    for i in range(coord,coord+200) :
      if i > chrlen[chrid] :
        break
      if i in picos['-'] :
        distancia = i-coord 
        suma += distancia
        count += 1
        break
 
mean_sep = suma/float(count)
offset = int(mean_sep/2+0.5) 
print("\tShifting value %d (from %d measures)" % (offset,count))  
print("\tTime taken: %d seconds" % (time.time()-time0))   

time0 = int(time.time())
print("Signal integration from strand signals")
depth_c = {}  
for chrid in chrseq.keys() :
  depth_c[chrid] = []
  depth[chrid]['+'] = [0]*offset + depth[chrid]['+'] 
  depth[chrid]['+'][-offset:] = [] 
  depth[chrid]['-'] = depth[chrid]['-']+[0]*offset  
  depth[chrid]['-'][:offset] = [] 
  for i in range(len(depth[chrid]['+'])) :
    depth_c[chrid].append(depth[chrid]['+'][i]+depth[chrid]['-'][i])
print("\tTime taken: %d seconds" % (time.time()-time0))   
#
if wig_int:
  time0 = int(time.time())
  print("Writing intermediate files: Integrated depth coverage")
  suffix = 'depth_c'
  WIG = open(WIG_file + "_%s.wig" % suffix,"w")
  for chrid in depth_c.keys() :
    WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    for s in depth_c[chrid] :
      WIG.write('%.4f\n' % s)
  WIG.close()
  print("\tTime taken: %d seconds" % (time.time()-time0))   

#
time0 = int(time.time())
print("Depth coverage wavelet denoising")
depth_c_wl = {} 
for chrid in depth_c.keys() :
  depth_c_wl[chrid] = denoising(depth_c[chrid])
print("\tTime taken: %d seconds" % (time.time()-time0))

time0 = int(time.time())
print("Mean depth coverage calculation")
if wig_int:
  print("Writing intermediate files: depth coverage wavelet denoised")
suma = 0
count = 0
suffix = 'depth_c_wl'
if wig_int: WIG = open(WIG_file + "_%s.wig" % suffix,"w")
for chrid in depth_c.keys() :
  if wig_int: WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
  if wig_int: WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
  for i in range(len(depth_c_wl[chrid])) :
    if depth_c_wl[chrid][i] <= 0 :
      depth_c_wl[chrid][i] = 0.0 
    else : 
      count += 1
      suma += depth_c_wl[chrid][i] 
    if wig_int: WIG.write('%.4f\n' % depth_c_wl[chrid][i])
if wig_int: WIG.close()
mean_c_wl = suma/count 
print("\tMean depth coverage: %.0f" % mean_c_wl)
print("\tTime taken: %d seconds" % (time.time()-time0))

#
time0 = int(time.time())
print("Writing results file: depth coverage wavelet denoised and normalized")
suffix = 'depth_c_wl_norm'
WIG = open(WIG_file + "_%s.wig" % suffix,"w")
for chrid in depth_c.keys() :
  WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
  WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
  for s in depth_c_wl[chrid] :
    WIG.write('%.4f\n' % (s/mean_c_wl))
WIG.close()
print("\tTime taken: %d seconds" % (time.time()-time0))

