library(psych)
library(tidyverse)
library(forcats)
library(viridis)

main = function(in.table, upstream, dnstream, cutoffpct, trimpct, ylab, refptlab, heatmap.cmap,
                meta.pal, avghmap.cmap, out.hmapsample, out.hmapgroup, out.metasample, out.metaoverlaysample,
                out.metaoverlayallsample, out.metagroup, out.metaoverlaygroup, out.metahmapsample, out.metahmapgroup){
    raw = read_tsv(in.table,
                   col_names=c("group","sample","index","position","cpm"),
                   col_types=cols(group=col_character(), sample=col_character(),
                                  index=col_integer(), position=col_double(),
                                  cpm=col_double()))
    raw$sample = fct_inorder(raw$sample, ordered = TRUE)
    raw$group = fct_inorder(raw$group, ordered = TRUE)
    
    nindices =  max(raw$index, na.rm=TRUE)
    nsamples = length(fct_unique(raw$sample))
    ngroups = length(fct_unique(raw$group))
    
    #percentile cutoff for heatmap visualization
    cutoff = quantile(raw$cpm, probs=cutoffpct, na.rm=TRUE)
    
    #plot heatmap facetted by sample and group
    heatmap_base = ggplot(data = raw %>% mutate_at(vars(cpm), funs(pmin(cutoff, .)))) +
        geom_raster(aes(x=position, y=index, fill=cpm)) +
        scale_y_reverse(name=paste(nindices, ylab), expand=c(0.01,0)) +
        scale_x_continuous(breaks = c(-upstream/1000, 0, dnstream/1000),
                           labels=c(ifelse(upstream>200, -upstream/1000, ''),
                                    refptlab,
                                    ifelse(dnstream>200, dnstream/1000, '')),
                           minor_breaks = scales::pretty_breaks(n=10),
                           name=paste("distance from", refptlab, "(kb)")) +
        scale_fill_viridis(option = heatmap.cmap, na.value="#FFFFFF00",
                           name="MNase-seq signal",
                           guide=guide_colorbar(title.position="top", barwidth=15,
                                                barheight=1, title.hjust=0.5)) +
        theme_minimal() +
        theme(text = element_text(size=12, face="bold", color="black"),
              legend.position = "top",
              legend.text = element_text(size=8, face="plain"),
              strip.text = element_text(size=12, face="bold", color="black"),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=12, face="bold", color="black", margin = unit(c(0,0,0,0),"cm")),
              panel.grid.major.x = element_line(color="black"),
              panel.grid.minor.x = element_line(color="grey80"),
              panel.grid.major.y = element_line(color="grey80"),
              panel.grid.minor.y = element_blank(),
              panel.spacing.x = unit(.5, "cm"))
  
    hmap.width = max(12, ((upstream+dnstream)/200)*(nsamples/ngroups))
    
    heatmap_samples = heatmap_base + facet_wrap(~sample, ncol=(nsamples/ngroups), dir="v")
    ggsave(out.hmapsample, plot = heatmap_samples,
           height= (.0005*nindices+7.5)*ngroups,
           width = hmap.width, units = "cm", limitsize=FALSE)
    rm(heatmap_samples)
    gc()
    heatmap_groups = heatmap_base + facet_wrap(~group, ncol=ngroups)
    ggsave(out.hmapgroup, plot = heatmap_groups,
           height= .002*nindices+14.75,
           width = hmap.width, units = "cm", limitsize=FALSE)
    rm(heatmap_groups)
    rm(heatmap_base)
    gc()
    
    #plot metagene and average heatmaps
    metadf.sample = raw %>% group_by(sample, position) %>%
        summarise(group = unique(group), mean = mean(cpm), sd = sd(cpm), trim.mean = mean(cpm, trim=trimpct),
                  win.mean = winsor.mean(cpm, trim=trimpct), win.sd = winsor.sd(cpm, trim=trimpct))
    
    meta.width = max(10, 2+(upstream+dnstream)/200)
    
    metaplot.sample= ggplot(data = metadf.sample,
                            aes(x=position, y=win.mean,
                                ymin=win.mean-win.sd, ymax=win.mean+win.sd,
                                color=group, fill=group)) +
                        geom_vline(xintercept=0, size=1) +
                        geom_ribbon(alpha=0.2, size=0) +
                        geom_line(size=1) +
                        scale_y_continuous(position="right", name=NULL) +
                        scale_x_continuous(name=paste("position relative to",
                                                      refptlab, "(kb)")) +
                        scale_fill_brewer(palette = meta.pal, guide=FALSE, direction=-1) +
                        scale_color_brewer(palette = meta.pal, guide=FALSE, direction=-1) +
                        facet_grid(sample~., switch="y") +
                        theme_minimal() +
                        theme(text = element_text(size=12, face="bold", color="black"),
                              axis.text.y = element_text(size=10, face="plain", color="black"),
                              axis.text.x = element_text(size=12, face="bold", color="black"),
                              strip.text.y = element_text(size=12, face="bold", color="black", angle=180, hjust=1),
                              panel.grid.major = element_line(color="grey85"),
                              panel.grid.minor = element_line(color="grey95"))
    
    ggsave(out.metasample, plot = metaplot.sample, height = 2*nsamples, width = meta.width, units = "cm")
    rm(metaplot.sample)
    gc()
    
    metaoverlay.sample= ggplot(data = metadf.sample,
                            aes(x=position, y=win.mean,
                                ymin=win.mean-win.sd, ymax=win.mean+win.sd,
                                group=sample, color=group, fill=group)) +
                        geom_vline(xintercept=0, size=1) +
                        geom_ribbon(alpha=0.05, size=0) +
                        geom_line(size=1, alpha=0.5) +
                        scale_y_continuous(position="right", name=NULL) +
                        scale_x_continuous(name=paste("position relative to",
                                                      refptlab, "(kb)")) +
                        scale_fill_brewer(palette = meta.pal, guide=FALSE, direction=-1) +
                        scale_color_brewer(palette = meta.pal, guide=FALSE, direction=-1) +
                        theme_minimal() +
                        facet_grid(group~., switch="y") +
                        theme(text = element_text(size=12, face="bold", color="black"),
                              axis.text.y = element_text(size=10, face="plain", color="black"),
                              axis.text.x = element_text(size=12, face="bold", color="black"),
                              strip.text.y = element_text(size=12, face="bold", color="black", angle=180, hjust=1),
                              panel.grid.major = element_line(color="grey85"),
                              panel.grid.minor = element_line(color="grey95"))
    
    ggsave(out.metaoverlaysample, plot = metaoverlay.sample, height = 3*ngroups, width = meta.width, units = "cm")
    
    metaoverlay.allsample= ggplot(data = metadf.sample,
                            aes(x=position, y=win.mean,
                                group=sample, color=group, fill=group)) +
                        geom_vline(xintercept=0, size=1) +
                        geom_line(size=1, alpha=0.5) +
                        scale_y_continuous(name="MNase-seq signal") +
                        scale_x_continuous(name=paste("position relative to",
                                                      refptlab, "(kb)")) +
                        scale_fill_brewer(palette = meta.pal, guide=FALSE, direction=-1) +
                        scale_color_brewer(palette = meta.pal, guide=guide_legend(), direction=-1) +
                        theme_minimal() +
                        theme(text = element_text(size=12, face="bold", color="black"),
                              axis.text.y = element_text(size=10, face="plain", color="black"),
                              axis.text.x = element_text(size=12, face="bold", color="black"),
                              panel.grid.major = element_line(color="grey85"),
                              panel.grid.minor = element_line(color="grey95"),
                              legend.title = element_blank(),
                              legend.text = element_text(size=12, face="bold", color="black"),
                              legend.position="top")
    ggsave(out.metaoverlayallsample, plot = metaoverlay.allsample, height=7, width = meta.width, units = "cm")
    
    metadf.group = raw %>% group_by(group, position) %>% 
        summarise(mean = mean(cpm), sd = sd(cpm), trim.mean = mean(cpm, trim=trimpct),
                  win.mean = winsor.mean(cpm, trim=trimpct), win.sd = winsor.sd(cpm, trim=trimpct))
    
    metaplot.group = ggplot(data = metadf.group,
                            aes(x=position, y=win.mean,
                                ymin=win.mean-win.sd, ymax=win.mean+win.sd,
                                color=group, fill=group)) +
                        geom_vline(xintercept=0, size=1) +
                        geom_ribbon(alpha=0.2, size=0) +
                        geom_line(size=1) +
                        scale_y_continuous(position="right", name=NULL,
                                           breaks = scales::pretty_breaks(n=3)) +
                        scale_x_continuous(name=paste("position relative to",
                                                      refptlab, "(kb)")) +
                        scale_fill_brewer(palette = meta.pal, guide=FALSE, direction=-1) +
                        scale_color_brewer(palette = meta.pal, guide=FALSE, direction=-1) +
                        facet_grid(group~., switch="y") +
                        theme_minimal() +
                        theme(text = element_text(size=12, face="bold", color="black"),
                              axis.text.y = element_text(size=10, face="plain", color="black"),
                              axis.text.x = element_text(size=12, face="bold", color="black"),
                              strip.text.y = element_text(angle=180, hjust=1),
                              panel.grid.major = element_line(color="grey85"),
                              panel.grid.minor = element_line(color="grey95"))
    
    ggsave(out.metagroup, plot = metaplot.group, height = 3*ngroups, width = meta.width, units = "cm")
    rm(metaplot.group)
    gc()
    metaoverlay.group = ggplot(data = metadf.group,
                            aes(x=position, y=win.mean,
                                ymin=win.mean-win.sd, ymax=win.mean+win.sd,
                                color=group, fill=group)) +
                        geom_vline(xintercept=0, size=1) +
                        geom_ribbon(alpha=0.1, size=0) +
                        geom_line(alpha=0.8, size=1) +
                        scale_y_continuous(name="MNase-seq signal") +
                        scale_x_continuous(name=paste("position relative to",
                                                      refptlab, "(kb)")) +
                        scale_fill_brewer(palette = meta.pal, direction=-1) +
                        scale_color_brewer(palette = meta.pal, direction=-1) +
                        theme_minimal() +
                        theme(text = element_text(size=12, face="bold", color="black"),
                              axis.text.y = element_text(size=10, face="plain", color="black"),
                              axis.text.x = element_text(size=12, face="bold", color="black"),
                              strip.text.y = element_text(angle=180, hjust=1),
                              panel.grid.major = element_line(color="grey85"),
                              panel.grid.minor = element_line(color="grey95"),
                              legend.title = element_blank(),
                              legend.position = "top",
                              legend.text = element_text(size=12, face="bold", color="black"))
    
    ggsave(out.metaoverlaygroup, plot = metaoverlay.group, height = 7, width = meta.width, units = "cm")
    rm(metaoverlay.group)
    gc()
    
    #average heatmap by sample
    metahmap.sample = ggplot(data = metadf.sample, aes(x=position, y=fct_rev(sample), fill=win.mean)) +
                            geom_raster() +
                            scale_fill_viridis(option = avghmap.cmap, name="MNase-seq signal", guide=
                                               guide_colorbar(title.position="top", title.hjust=0.5, barwidth=8)) +
                            scale_x_continuous(expand = c(0.01, 0),
                                               name=paste("position relative to", refptlab, "(kb)")) +
                            theme_minimal() +
                            theme(text = element_text(size=12, face="bold", color="black"),
                                  legend.position="top",
                                  axis.title.y = element_blank(),
                                  axis.text = element_text(size=12, face="bold", color="black"),
                                  legend.text = element_text(size=10, face="plain"),
                                  panel.grid.major.x = element_line(color="black", size=1),
                                  panel.grid.major.y = element_blank())
    ggsave(out.metahmapsample, plot = metahmap.sample, height = nsamples+2, width = meta.width, units = "cm")
    rm(metahmap.sample)
    gc()
    
    #average heatmap by group
    metahmap.group = ggplot(data = metadf.group, aes(x=position, y=fct_rev(group), fill=win.mean)) +
                            geom_raster() +
                            scale_fill_viridis(option = avghmap.cmap, name="MNase-seq signal", guide=
                                               guide_colorbar(title.position="top", title.hjust=0.5, barwidth=8)) +
                            scale_x_continuous(expand = c(0.01, 0),
                                               name=paste("position relative to", refptlab, "(kb)")) +
                            theme_minimal() +
                            theme(text = element_text(size=12, face="bold", color="black"),
                                  legend.position="top",
                                  axis.title.y = element_blank(),
                                  axis.text = element_text(size=12, face="bold", color="black"),
                                  legend.text = element_text(size=10, face="plain"),
                                  panel.grid.major.x = element_line(color="black", size=1),
                                  panel.grid.major.y = element_blank())
    
    ggsave(out.metahmapgroup, plot = metahmap.group, height = 1.5*ngroups+2, width = meta.width, units = "cm")
    rm(metahmap.group)
    gc()
}

main(
    in.table = snakemake@input[["matrix"]],
    upstream = snakemake@params[["upstream"]],
    dnstream = snakemake@params[["dnstream"]],
    cutoffpct = snakemake@params[["pct_cutoff"]],
    trimpct = snakemake@params[["trim_pct"]],
    ylab = snakemake@params[["ylabel"]], 
    refptlab = snakemake@params[["refpointlabel"]],
    heatmap.cmap = snakemake@params[["heatmap_cmap"]],
    meta.pal = snakemake@params[["metagene_palette"]],
    avghmap.cmap = snakemake@params[["avg_heatmap_cmap"]],
    out.hmapsample = snakemake@output[["heatmap_sample"]],
    out.hmapgroup = snakemake@output[["heatmap_group"]],
    out.metasample = snakemake@output[["metagene_sample"]],
    out.metaoverlaysample = snakemake@output[["metagene_sample_overlay"]],
    out.metaoverlayallsample = snakemake@output[["metagene_sample_overlay_all"]],
    out.metagroup = snakemake@output[["metagene_group"]],
    out.metaoverlaygroup = snakemake@output[["metagene_overlay"]],
    out.metahmapsample = snakemake@output[["metaheatmap_sample"]],
    out.metahmapgroup = snakemake@output[["metaheatmap_group"]]
)
