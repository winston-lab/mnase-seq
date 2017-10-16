library(tidyverse)
library(forcats)
library(ggridges)

raw = read_tsv(snakemake@input[["table"]], col_names = c("sample", "group", "fragsize"))
raw$sample = as_factor(raw$sample) %>% fct_rev()

#p = ggplot(data = raw, aes(x=sample, y=fragsize, fill=group)) +
#        geom_violin(scale="count", adjust=1, draw_quantiles=c(0.5)) +
#        scale_fill_brewer(palette = "Set1") + 
#        ylab("fragment size (bp)") +
#        coord_flip() +
#        theme_bw() +
#        theme(axis.text = element_text(size=14),
#              axis.title = element_text(size=14, face="bold"))

p2 = ggplot(data = raw, aes(x=fragsize, y=sample, fill=group)) +
        geom_density_ridges() +
        scale_fill_brewer(palette="Set1", guide=FALSE) +
        scale_x_continuous(expand = c(0.01, 0),
                           name="fragment size (bp)",
                           breaks = scales::pretty_breaks(n=10)) +
        scale_y_discrete(expand = c(0.01, 0)) +
        theme_minimal() +
        theme(text = element_text(size=12, face="bold", color="black"),
              axis.text = element_text(size=12, face="bold", color="black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(vjust=0))

h = length(fct_unique(raw$sample))

ggsave(snakemake@output[["plot"]], plot=p2, width = 10, height=h, units="cm")
