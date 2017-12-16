library(tidyverse)

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

melt(inmatrix = snakemake@input[["matrix"]],
     refpt = snakemake@params[["refpoint"]],
     group = snakemake@params[["group"]],
     sample = snakemake@wildcards[["sample"]],
     binsize = snakemake@params[["binsize"]],
     upstream = snakemake@params[["upstream"]],
     outpath = snakemake@output[[1]])
