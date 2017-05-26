library(tidyverse)
library(forcats)

raw = read_table2(gzfile(snakemake@input[[1]]), col_names = c("sample", "count"))

raw$sample = as.factor(raw$sample) %>% fct_rev()

lim = quantile(raw$count, .995)

violin = ggplot(data = raw, aes(x=sample, y=count)) +
          geom_violin(bw=.6, draw_quantiles = .5) +
          scale_y_continuous(limits = c(0, lim)) +
          ylab("depth (fragments per base)") +
          coord_flip() +
          theme_bw() +
          theme(axis.title = element_text(size=14, face="bold"), 
                axis.text = element_text(size=14))

h = 3*length(fct_unique(raw$sample))

ggsave(snakemake@output[[1]], plot = violin, width = 10, height=h, units="cm")
          
