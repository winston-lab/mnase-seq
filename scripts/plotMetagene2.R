library(tidyverse)
library(viridis)

raw = read_table2(snakemake@input[[1]], col_names=c("sample", "position","cpm") ) %>% filter(position <=1)

knots = round((max(raw$position) - min(raw$position))*1000/148*3)

p = ggplot(data = raw, aes(x=position, y=cpm, color=sample, fill=sample)) +
    geom_vline(xintercept = 0, size=2) +
    geom_smooth(method="gam", formula = y ~ s(x, bs="cr", k=knots), size=1.5, alpha=0.6) +
    xlab(paste("position relative to", snakemake@params[["refpointlabel"]], "(kb)")) +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme_minimal() +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12, face="bold"))

w = round((max(raw$position) - min(raw$position))*1000/148)*1.5

ggsave(snakemake@output[[1]], plot = p, width = w, height = 10, units = "cm")
