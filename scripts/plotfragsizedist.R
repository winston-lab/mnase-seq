library(tidyverse)
library(forcats)

raw = read_tsv(snakemake@input[["table"]], col_names = c("sample", "fragsize"))
raw$sample = as.factor(raw$sample) %>% fct_rev()

p = ggplot(data = raw, aes(x=sample, y=fragsize)) +
      geom_violin(scale="count", adjust=1, draw_quantiles=c(0.5)) +
      ylab("fragment size (bp)") +
      coord_flip() +
      theme_bw() +
      theme(axis.text = element_text(size=14),
            axis.title = element_text(size=14, face="bold"))

h = 3*length(fct_unique(raw$sample))

ggsave(snakemake@output[["plot"]], plot=p, width = 10, height=h, units="cm")
