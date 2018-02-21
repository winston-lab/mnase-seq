#!/usr/bin/env/ Rscript
library(argparse)
library(tidyverse)

parser = ArgumentParser()
parser$add_argument('-i', '--input', type='character')
parser$add_argument('-r', '--refpt', type='character', nargs="+")
parser$add_argument('-g', '--group', type='character', nargs="+")
parser$add_argument('-s', '--sample', type='character', nargs="+")
parser$add_argument('-a', '--annotation', type='character', nargs="+")
parser$add_argument('-b', '--binsize', type='integer')
parser$add_argument('-u', '--upstream', type='integer')
parser$add_argument('-o', '--output', type='character')

args = parser$parse_args()

melt = function(inmatrix, refpt, group, sample, annotation, binsize, upstream, outpath){
    raw = read_tsv(inmatrix, skip=3, col_names=FALSE)
    names(raw) = seq(ncol(raw))
    
    df = raw %>%
          rownames_to_column(var="index") %>%
          gather(key = variable, value=value, -index, convert=TRUE) %>%
          filter(!is.na(value)) 
    if(binsize>1){
        df = df %>%
            transmute(group = group, sample = sample,
                      annotation = annotation, 
                      index = as.integer(index),
                      position = (as.numeric(variable)*binsize-(upstream+1.5*binsize))/1000,
                      cpm = as.numeric(value)) 
    } else if (refpt=="TES"){
        df = df %>%
            transmute(group = group, sample = sample,
                      annotation = annotation, 
                      index = as.integer(index),
                      position = (as.numeric(variable)-(1+upstream))/1000,
                      cpm = as.numeric(value)) 
    } else {
        df = df %>%
            transmute(group = group, sample = sample,
                      annotation = annotation, 
                      index = as.integer(index),
                      position = (as.numeric(variable)-(2+upstream))/1000,
                      cpm = as.numeric(value)) 
    }
    write_tsv(df, path=outpath, col_names=FALSE)
    return(df)
}

melt(inmatrix = args$input,
     refpt = args$refpt,
     group = paste(args$group, collapse=" "),
     sample = paste(args$sample, collapse=" "),
     annotation = paste(args$annotation, collapse=" "),
     binsize = args$binsize,
     upstream = args$upstream,
     outpath = args$output)
