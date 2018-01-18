#!/usr/bin/env/ Rscript
library(optparse)
library(tidyverse)

option_list = list(make_option(c("-i", "--input")),
                   make_option(c("-r", "--refpt")),
                   make_option(c("-g", "--group")),
                   make_option(c("-s", "--sample")),
                   make_option(c("-b", "--binsize"), type="integer"),
                   make_option(c("-u", "--upstream"), type="integer"),
                   make_option(c("-o", "--output")))

opt = parse_args(OptionParser(option_list=option_list))

melt = function(inmatrix, refpt, group, sample, binsize, upstream, outpath){
    raw = read_tsv(inmatrix, skip=3, col_names=FALSE)
    names(raw) = seq(ncol(raw))
    
    df = raw %>%
          rownames_to_column(var="index") %>%
          gather(key = variable, value=value, -index, convert=TRUE) %>%
          filter(!is.na(value)) 
    if(binsize>1){
        df = df %>%
            transmute(group = group,
                      sample = sample,
                      index = as.integer(index),
                      position = (as.numeric(variable)*binsize-(upstream+1.5*binsize))/1000,
                      cpm = as.numeric(value)) 
    } else if (refpt=="TES"){
        df = df %>%
            transmute(group = group,
                      sample = sample,
                      index = as.integer(index),
                      position = (as.numeric(variable)-(1+upstream))/1000,
                      cpm = as.numeric(value)) 
    } else {
        df = df %>%
            transmute(group = group,
                      sample = sample,
                      index = as.integer(index),
                      position = (as.numeric(variable)-(2+upstream))/1000,
                      cpm = as.numeric(value)) 
    }
    write_tsv(df, path=outpath, col_names=FALSE)
    return(df)
}

melt(inmatrix = opt$input,
     refpt = opt$refpt,
     group = opt$group,
     sample = opt$sample,
     binsize = opt$binsize,
     upstream = opt$upstream,
     outpath = opt$output)
