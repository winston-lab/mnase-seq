library(tidyverse)
library(reshape2)

raw = read_table2(snakemake@input[["matrix"]], skip=3, col_names=FALSE)

binsize = snakemake@params[["binsize"]]
upstream = snakemake@params[["upstream"]]
downstream = snakemake@params[["dnstream"]]
n = nrow(raw)

df = raw %>%
      rownames_to_column(var="index") %>%
      melt(id.vars="index") %>%
      as_tibble() %>%
      transmute(sample = snakemake@params[["name"]], position = (as.numeric(variable)*binsize-upstream)/1000, cpm = as.numeric(value))

write.table(df, file=gzfile(snakemake@output[[1]]), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
