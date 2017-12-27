library(tidyverse)
library(forcats)
library(viridis)
library(dendsort)

format_xaxis_kb = function(refptlabel){
    function(x) if_else(x==0, refptlabel, as.character(x))
}

format_xaxis_nt = function(refptlabel){
    function(x) if_else(x==0, refptlabel, as.character(x*1000))
}

label_xaxis= function(ggp, refptlabel, upstream, dnstream){
    if(upstream>1000 | dnstream>1000){
        ggp = ggp +
            scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                               labels=format_xaxis_kb(refptlabel=refptlabel),
                               name=paste("distance from", refptlabel, "(kb)"),
                               limits = c(-upstream/1000, dnstream/1000),
                               expand=c(0.05,0))
    } else{
        ggp = ggp +
            scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                               labels=format_xaxis_nt(refptlabel=refptlabel),
                               name=paste("distance from", refptlabel, "(nt)"),
                               limits = c(-upstream/1000, dnstream/1000),
                               expand=c(0.05,0))
    } 
    return(ggp)
}

hmap = function(df, nindices, ylabel, upstream, dnstream, refptlab, cmap){
    heatmap_base = ggplot(data = df) +
        geom_raster(aes(x=position, y=index, fill=cpm)) +
        scale_y_reverse(name=paste(nindices, ylabel), expand=c(0.02, 0)) +
        scale_fill_viridis(option = cmap, na.value="FFFFFF00", name='MNase-seq signal', guide=guide_colorbar(title.position="top", barwidth=15, barheight=1, title.hjust=0.5)) +
        theme_minimal() +
        theme(text = element_text(size=12, face="bold", color="black"),
              legend.position = "top",
              legend.title = element_text(size=12, face="bold", color="black"),
              legend.text = element_text(size=10, face="plain"),
              legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(0,0,0,0),
              strip.text = element_text(size=12, face="bold", color="black"),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=12, face="bold", color="black", margin = unit(c(0,0,0,0),"cm")),
              panel.grid.major.x = element_line(color="black"),
              panel.grid.minor.x = element_line(color="black"),
              panel.grid.major.y = element_line(color="black"),
              panel.grid.minor.y = element_blank(),
              panel.spacing.x = unit(.5, "cm"))
    heatmap_base = heatmap_base %>%
        label_xaxis(refptlabel=refptlab, upstream=upstream, dnstream=dnstream)
    return(heatmap_base)
}
    
main = function(intable, samplelist, upstream, dnstream, pct_cutoff,
                cluster, k, refptlab, ylabel, cmap, samples_out, group_out){
    raw = read_tsv(intable, col_names=c("group", "sample", "index", "position","cpm")) %>%
        filter(sample %in% samplelist & !is.na(cpm)) %>% 
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))
        
    nindices = max(raw$index, na.rm=TRUE)
    nsamples = length(samplelist)
    ngroups = length(fct_unique(raw$group))
    
    #clustering
    if (cluster){
        #first k-means clustering on NOTE: ??unscaled?? data
        rr = raw %>% select(-group) %>% unite(cid, c(sample, position), sep="~") %>%
                spread(cid, cpm, fill=0) %>% select(-index)
        clust = kmeans(rr, centers=k)
       
        #then hierarchical clustering on the k-means centers and sorting
        #of the resulting dendrogram
        centerclust = clust$centers %>% dist() %>% hclust() %>% dendsort(isReverse=TRUE)
        
        reorder = clust$cluster %>% as_tibble() %>% rename(cluster=value) %>%
            mutate_at(vars(cluster), funs(factor(., levels=centerclust$order, ordered=TRUE))) %>%
            mutate(og_index=row_number()) %>%
            arrange(cluster,og_index) %>%
            mutate(new_index=row_number())
        raw = raw %>% left_join(reorder, by=c("index"="og_index")) %>%
            select(-index) %>% rename(index=new_index)
        rm(rr, clust, centerclust, reorder)
    } 
    
    #get replicate info for sample facetting
    repl_df = raw %>% select(group, sample) %>% distinct() %>% group_by(group) %>%
        mutate(replicate=row_number()) %>% ungroup() %>% select(-group)
    
    #plot heatmap facetted by sample and group
    df_sample = raw %>% left_join(repl_df, by="sample")
    sample_cutoff = quantile(df_sample$cpm, probs=pct_cutoff, na.rm=TRUE)
    df_sample = df_sample %>%
        mutate_at(vars(cpm),
                  funs(pmin(sample_cutoff, .)))
    hmap_sample = hmap(df=df_sample, nindices=nindices, ylabel=ylabel,
                       upstream=upstream, dnstream=dnstream, refptlab=refptlab,
                       cmap=cmap) +
        facet_grid(replicate~group) +
        theme(strip.text.y=element_text(angle=0))
 
    hmap.width = max(12, (.0008*(upstream+dnstream)+3.4)*ngroups) 
    ggsave(samples_out, plot=hmap_sample, height=(.0005*nindices+7.5)*max(repl_df$replicate),
           width=hmap.width, units="cm", limitsize=FALSE)
    rm(hmap_sample, df_sample, repl_df)

    
    df_group = raw %>% group_by(group, index, position) %>% summarise(cpm = mean(cpm))
    rm(raw)
    group_cutoff = quantile(df_group$cpm, probs=pct_cutoff, na.rm=TRUE)
    df_group = df_group %>% 
        mutate_at(vars(cpm),
                  funs(pmin(group_cutoff, .)))
    hmap_group = hmap(df=df_group, nindices=nindices, ylabel=ylabel,
                       upstream=upstream, dnstream=dnstream, refptlab=refptlab,
                       cmap=cmap) +
        facet_wrap(~group, ncol=ngroups)
    ggsave(group_out, plot = hmap_group, height= .0009*nindices+11.5,
           width=hmap.width, units="cm")
}

main(intable= snakemake@input[["matrix"]],
     samplelist = snakemake@params[["samplelist"]],
     upstream = snakemake@params[["upstream"]],
     dnstream= snakemake@params[["dnstream"]],
     pct_cutoff = snakemake@params[["pct_cutoff"]],
     cluster = snakemake@params[["cluster"]],
     k = snakemake@params[["nclust"]],
     refptlab = snakemake@params[["refpointlabel"]],
     ylabel = snakemake@params[["ylabel"]],
     cmap = snakemake@params[["heatmap_cmap"]],
     samples_out = snakemake@output[["heatmap_sample"]],
     group_out = snakemake@output[["heatmap_group"]])
