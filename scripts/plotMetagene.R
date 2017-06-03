library(tidyverse)
library(viridis)

raw = read_table2(snakemake@input[[1]], col_name=c("sample","position","cpm"))

p = ggplot(data = raw, aes(x=position, y=cpm, group=sample, color=sample, fill=sample)) +
      geom_vline(xintercept=0, size=2) +
      geom_smooth(method="gam", formula = y ~ s(x, bs="gp", k=30, m=c(1,.1)), size=1.5, alpha=.4) +
      #geom_point(alpha=0.1) +
      xlab(paste("distance to", snakemake@params[["refpointlabel"]], "(kb)")) +
      ylab("normalized midpoint occupancy") +
      #scale_y_continuous(limits = c(0, .6)) +
      scale_color_viridis(discrete=TRUE) +
      scale_fill_viridis(discrete=TRUE) +
      theme_minimal() +
      theme(legend.position = "right",
            #legend.position=c(.3,.92),
            legend.title = element_blank(),
            axis.text.x = element_text(size=12),
            axis.text.y = element_blank(),
            axis.title = element_text(size=12, face="bold"))

ggsave(snakemake@output[[1]], plot=p, width=14, height=8, units="cm")
