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






#p = ggplot(data = df) +
#  geom_tile(aes(x=position, y=index, fill=cpm)) +
#  #annotate("text", x=(-upstream/1000), y=1.02*n, label = paste((-upstream/1000), "\nkb", sep=""), size=3, hjust="center", vjust=1, fontface="bold") +
#  annotate("text", x=0, y=1.02*n, label = snakemake@params[["refpointlabel"]], size=4, fontface="bold", hjust="center", vjust=0) +
#  annotate("text", x=(downstream/1000), y=1.02*n, label = paste("+", (downstream/1000), "kb", sep=""), size=4, hjust="right", vjust=0, fontface="bold") +
#  scale_y_reverse(name=snakemake@params[["ylabel"]], breaks=NULL) +
#  scale_x_continuous(name=NULL, breaks = c(-upstream/1000, 0, downstream/1000), labels=NULL) +
#  scale_fill_viridis(na.value="white", option=snakemake@params[["cmap"]], name="normalized\nmidpoint\noccupancy") +
#  theme_void() 
  #theme(legend.position = c(.95,.9), legend.justification = c(1,1))

#ggsave(filename = snakemake@output[[1]], plot=p, width=snakemake@params[["figwidth"]], height=snakemake@params[["figheight"]], units = snakemake@params[["figunits"]])
