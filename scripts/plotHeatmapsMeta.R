library(tidyverse)
library(forcats)
library(viridis)
library(RColorBrewer)

raw = read_table2(snakemake@input[["matrix"]],
	 col_names=c("group", "sample", "index", "position","cpm"),
	 col_types=cols(group=col_character(), sample=col_character(), index=col_integer(), position=col_double(), cpm=col_double()))
raw$sample = factor(raw$sample, ordered = TRUE)
raw$group = factor(raw$group, ordered = TRUE)

#for cubic spline fitting, set number of knots to around 3 per expected nucleosome (i.e. per 148 bp)
knots = round((max(raw$position) - min(raw$position))*1000/148*3)
#number of positions to evaluate smoothed fit at
n_avg = 1000/snakemake@params[["binsize"]]*(max(raw$position)- min(raw$position))

nindices =  max(raw$index)
nsamples = length(fct_unique(raw$sample))
ngroups = length(fct_unique(raw$group))
w = round((max(raw$position) - min(raw$position))*1000/148)

upstream = snakemake@params[["upstream"]]
downstream = snakemake@params[["dnstream"]]

#pseudocount for log-transform
#pcount = .01

#percentile cutoff for heatmap visualization
cutoff = quantile(raw$cpm, probs=snakemake@params[["pct_cutoff"]], na.rm=TRUE)

#plot heatmap facetted by sample and group
heatmap_base = ggplot(data = raw %>% mutate_at(vars(cpm), funs(pmin(cutoff, .)))) + 
  geom_tile(aes(x=position, y=index, fill=cpm)) +
  scale_y_reverse(name=paste(nindices, snakemake@params[["ylabel"]])) +
  scale_x_continuous(breaks = c(-upstream/1000, 0, downstream/1000), labels=c(-upstream/1000, snakemake@params[["refpointlabel"]], downstream/1000)) +
  xlab(paste("distance from", snakemake@params[["refpointlabel"]], "(kb)")) +
  scale_fill_viridis(option = snakemake@params[["heatmap_cmap"]],na.value="white", name="normalized MNase-seq fragment midpoints", guide=guide_colorbar(title.position="top", barwidth=15, barheight=1, title.hjust=0.5)) +
  theme_minimal() +
  theme(strip.text = element_text(size=12, face="bold"),
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.ticks.length = unit(-2, "mm"))


heatmap_samples = heatmap_base + facet_wrap(~sample, ncol=(nsamples/ngroups))
ggsave(snakemake@output[["heatmap_sample"]], plot = heatmap_samples, height=10+round((nindices/600)*(nsamples/ngroups)), width = 8+.3*w*ngroups, units = "cm")
rm(heatmap_samples)
heatmap_groups = heatmap_base + facet_wrap(~group, ncol=(nsamples/ngroups))
ggsave(snakemake@output[["heatmap_group"]], plot = heatmap_groups, height=10+round(3*(nindices/1000)), width = 8+.3*w*ngroups, units = "cm")
rm(heatmap_groups)
rm(heatmap_base)

#plot metagene and average heatmaps
#metagene_base = ggplot(data = raw %>% filter(cpm <= cutoff), aes(x=position, y=cpm, color=group, fill=group)) +
metagene_base = ggplot(data = raw, aes(x=position, y=cpm, color=group, fill=group)) +
  geom_vline(xintercept = 0, size=2) +
  geom_smooth(method="gam", formula = y ~ s(x, bs="cr", k=knots), size=1.5, alpha=0.6, n=n_avg) +
  scale_y_continuous(position="right") +
  xlab(paste("position relative to", snakemake@params[["refpointlabel"]], "(kb)")) +
  scale_fill_brewer(palette=snakemake@params[["metagene_palette"]]) +
  scale_color_brewer(palette=snakemake@params[["metagene_palette"]]) +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        #axis.title.y.right = element_text(angle=0, vjust=0.5),
        axis.title.y.right = element_blank(),
        strip.text.y = element_text(size=12, angle=180, face="bold", hjust=1))

metagene_samples = metagene_base + facet_grid(sample~., switch="y") + theme(legend.position = "none")
ggsave(snakemake@output[["metagene_sample"]], plot = metagene_samples, height = 3*nsamples, width = 5+.7*w, units = "cm")

#average heatmap by sample
metagene_samples.fit = ggplot_build(metagene_samples)$data[[2]]
rm(metagene_samples)
metagene_samples.fit$sample = factor(metagene_samples.fit$PANEL, labels = levels(raw$sample), ordered = TRUE) %>% fct_rev

avg_heatmap_samples = ggplot(data = metagene_samples.fit, aes(x=x, y=sample)) +
  geom_tile(aes(fill=y)) +
  scale_fill_viridis(option = snakemake@params[["avg_heatmap_cmap"]], guide_colorbar(title = "normalized\nMNase-seq\nfragment\nmidpoints")) +
  xlab(paste("position relative to", snakemake@params[["refpointlabel"]], "(kb)")) +
  theme_minimal() +
  theme(axis.text = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face="bold"),
        legend.title = element_text())
ggsave(snakemake@output[["metaheatmap_sample"]], plot = avg_heatmap_samples, height = 2*nsamples, width = 8+.5*w, units = "cm")
rm(avg_heatmap_samples)

#metagene by group
metagene_group = metagene_base + facet_grid(group~., switch="y") + theme(legend.position = "none")
ggsave(snakemake@output[["metagene_group"]], plot = metagene_group, height = 4*ngroups, width = 5+.7*w, units = "cm")

#average heatmap by group
metagene_group.fit = ggplot_build(metagene_group)$data[[2]]
rm(metagene_group)
metagene_group.fit$group = factor(metagene_group.fit$group, labels = levels(raw$group), ordered = TRUE) %>% fct_rev

avg_heatmap_groups = ggplot(data = metagene_group.fit, aes(x=x, y=group)) +
  geom_tile(aes(fill=y)) +
  scale_fill_viridis(option = snakemake@params[["avg_heatmap_cmap"]], guide_colorbar(title = "normalized\nMNase-seq\nfragment\nmidpoints")) +
  xlab(paste("position relative to", snakemake@params[["refpointlabel"]], "(kb)")) +
  theme_minimal() +
  theme(axis.text = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face="bold"),
        legend.title = element_text())
ggsave(snakemake@output[["metaheatmap_group"]], plot = avg_heatmap_groups, height = 2*ngroups, width = 8+.5*w, units = "cm")
rm(avg_heatmap_groups)

#plot overlaid metagene
metagene_overlay = metagene_base + theme(legend.position = "left", legend.title = element_blank(), legend.text = element_text(face="bold"))
ggsave(snakemake@output[["metagene_overlay"]], plot = metagene_overlay, height = 10, width = 8+.8*w, units = "cm")
rm(metagene_overlay)
rm(metagene_base)
