library(tidyverse)

melt = function(inmatrix, group, sample, binsize, upstream, outpath){
    raw = read_tsv(inmatrix, skip=3, col_names=FALSE)
    names(raw) = seq(ncol(raw))
    
    df = raw %>%
          rownames_to_column(var="index") %>%
          gather(key = variable, value=value, -index, convert=TRUE) %>%
          filter(!is.na(value)) %>%
          transmute(group = group, 
                    sample = sample,
                    index = as.numeric(index),
                    position = (as.numeric(variable)*binsize-upstream)/1000,
                    cpm = as.numeric(value)) 
    
    write_tsv(df, path=outpath, col_names=FALSE)
    return(df)
}

melt(inmatrix = snakemake@input[["matrix"]],
     group = snakemake@params[["group"]],
     sample = snakemake@wildcards[["sample"]],
     binsize = snakemake@params[["binsize"]],
     upstream = snakemake@params[["upstream"]],
     outpath = snakemake@output[[1]])
