library(tidyverse)
library(forcats)
library(viridis)
library(ggthemes)
library(ggrepel)

import = function(path){
    read_tsv(path) %>%
        mutate_at(vars(sample, status), funs(fct_inorder(., ordered=TRUE))) %>%
        return()
}
            
main = function(seq_len_dist_in, per_tile_in, per_base_qual_in,
                per_base_seq_in, per_base_n_in, per_seq_gc_in,
                per_seq_qual_in, adapter_content_in, seq_dup_in, kmer_in,
                seq_len_dist_out, per_tile_out, per_base_qual_out,
                per_base_seq_out, per_seq_gc_out, per_seq_qual_out,
                adapter_content_out, seq_dup_out, kmer_out){
    #damnit fastqc...why bin some of the data and then output in this shite format...
    length_distribution = import(seq_len_dist_in) %>% 
        separate(length, into=c('a','b'), sep="-", fill="right", convert=TRUE) %>% 
        mutate_at(vars(count), funs(if_else(is.na(b), ., ./2))) %>% 
        gather(key, length, c(a,b)) %>% 
        filter(!is.na(length)) %>% 
        select(-key)
    
    nsamples = n_distinct(length_distribution$sample)
    
    per_tile_quality = import(per_tile_in) %>% 
        mutate_at(vars(tile), funs(fct_inorder(as.character(.), ordered=TRUE)))
    
    per_base_qual = import(per_base_qual_in) %>%
        left_join(length_distribution, by=c("base"="length", "sample", "status")) %>% 
        group_by(sample, status) %>%
        mutate_at(vars(count), funs(if_else(is.na(.), 0, .))) %>% 
        mutate(n = lag(sum(count)-cumsum(count), default=sum(count))) %>% 
        mutate_at(vars(n), funs(./max(n)))
    
    adapter_content = import(adapter_content_in)
    
    per_base_seq = import(per_base_seq_in) %>%
        left_join(import(per_base_n_in), by=c("base", "sample","status")) %>% 
        rename(position=base, n=n_count) %>% 
        gather(base, pct, -c(position, sample, status)) %>% 
        mutate_at(vars(base), funs(toupper(.))) %>% 
        mutate_at(vars(base), funs(fct_inorder(., ordered=TRUE)))
    
    per_seq_gc = import(per_seq_gc_in) %>%
        filter(count > 0) %>% 
        group_by(sample, status) %>% 
        mutate(norm_count = count/sum(count))
    
    per_seq_qual = import(per_seq_qual_in) %>% 
        filter(count > 0) %>% 
        group_by(sample, status) %>% 
        mutate(norm_count = count/sum(count))
    
    duplication_levels = import(seq_dup_in) %>% 
        mutate_at(vars(duplication_level), funs(fct_inorder(., ordered=TRUE)))
    
    kmer_content = import(kmer_in)
    
    theme_standard = theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black"),
              axis.title = element_text(size=12, color="black", face="bold"),
              strip.placement = "outside",
              strip.background = element_blank(),
              strip.text = element_text(size=12, color="black", face="bold"),
              strip.text.y = element_text(angle=-180, hjust=1))
    
    length_dist_plot = ggplot(data = length_distribution %>%
                                     group_by(sample, status) %>%
                                     mutate(normcount = count/max(count)),
                                 aes(x=length, y=normcount)) +
        geom_col(fill="#114477") +
        scale_x_continuous(breaks=scales::pretty_breaks(n=6), name="read length (nt)") +
        scale_y_continuous(breaks=scales::pretty_breaks(n=2), name="normalized counts") +
        facet_grid(sample~status, scales="free_y", switch="y") +
        ggtitle("read length distributions") +
        theme_standard +
        theme(axis.text.y = element_text(size=10, face="plain"))
    
    ggsave(seq_len_dist_out, plot=length_dist_plot, width=26, height=2+2*nsamples, units="cm")
    
    tile_quality_plot = ggplot(data = per_tile_quality %>% filter(status=="raw"),
                               aes(x=base, y=tile, fill=mean)) +
        geom_raster() +
        scale_fill_viridis(direction=-1, guide=guide_colorbar(title="mean\nquality\nscore", barheight=10)) +
        scale_x_continuous(expand=c(0,0), name="cycle number", breaks=scales::pretty_breaks(n=6)) +
        ylab("flow cell tile") +
        ggtitle("per tile sequencing quality") +
        facet_grid(sample~status, scales="free_y", switch="y") +
        theme_standard +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major.y = element_blank(),
              strip.text.x = element_blank())
    
    ggsave(per_tile_out, plot=tile_quality_plot, width=14, height=2+2*nsamples, units="cm")
    
    per_base_qual_plot = ggplot(data=per_base_qual, aes(x=base, y=fct_rev(sample),
                                               height=n, fill=mean, color=mean)) +
        geom_tile() +
        scale_color_viridis(guide=FALSE, direction=-1) +
        scale_fill_viridis(guide=guide_colorbar(title="mean quality", barwidth=12,
                                                barheight=1, title.position = "top",
                                                title.hjust=0.5),
                           direction=-1) +
        scale_x_continuous(expand=c(0,0), name="read length (nt)", breaks=scales::pretty_breaks(n=6)) +
        facet_grid(sample~status, scales="free_y", switch="y") +
        ggtitle("per base sequencing quality",
                subtitle = expression("bar height " %prop% " fraction of reads")) +
        theme_standard +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              plot.subtitle = element_text(size=12),
              legend.position="top",
              legend.margin = margin(0,0,0,0))
    
    ggsave(per_base_qual_out, plot=per_base_qual_plot, width=26, height=2.5+1.25*nsamples, units="cm")
    
    adapter_plot = ggplot(data = adapter_content, aes(x=position, y=0, fill=pct, color=pct)) +
        geom_raster() +
        scale_color_viridis(guide=FALSE) +
        scale_fill_viridis(guide=guide_colorbar(title="% reads with adapter", barwidth=12,
                                                barheight=1, title.position = "top",
                                                title.hjust=0.5)) +
        scale_x_continuous(expand=c(0,0), name="read length (nt)") +
        scale_y_continuous(expand=c(0,0), name=NULL, breaks=0, labels=NULL) +
        facet_grid(sample~status, switch="y") +
        ggtitle("adapter content") +
        theme_standard +
        theme(legend.position="top",
              legend.margin = margin(0,0,0,0))
    
    ggsave(adapter_content_out, plot=adapter_plot, width=26, height=2+1.25*nsamples, units="cm")
              
    per_base_seq_plot = ggplot(data = per_base_seq, aes(x=position, y=pct, color=base)) +
        geom_line() +
        scale_color_ptol(guide=guide_legend(label.position="top", label.hjust=0.5,
                                            keyheight=0.2)) +
        scale_x_continuous(expand=c(0,0), name="position in read", breaks=scales::pretty_breaks(n=6)) +
        scale_y_continuous(name="% of reads", breaks=scales::pretty_breaks(n=2)) +
        facet_grid(sample~status, switch="y") +
        ggtitle("per base sequence content") +
        theme_standard +
        theme(legend.position="top",
              legend.title = element_blank(),
              legend.margin = margin(0,0,0,0),
              legend.key.size = unit(1, "cm"),
              legend.text = element_text(size=12, face="bold"),
              axis.text.y = element_text(size=10, face="plain"))
    
    ggsave(per_base_seq_out, plot=per_base_seq_plot, width=26, height=2+2.25*nsamples, units="cm")
    
    per_seq_gc_plot = ggplot(data = per_seq_gc, aes(x=gc_content, y=norm_count)) +
        geom_line(color="#114477") +
        scale_x_continuous(expand=c(0,0), name="GC%") +
        #xlab("GC%") +
        scale_y_continuous(breaks=scales::pretty_breaks(n=2), name="normalized counts") +
        facet_grid(sample~status, switch="y") +
        ggtitle("per sequence GC content") +
        theme_standard +
        theme(axis.text.y = element_text(size=10, face="plain"),
              panel.spacing.x = unit(1, "cm"),
              plot.margin = margin(5.5, 12, 5.5, 5.5, unit="pt"))
    
    ggsave(per_seq_gc_out, plot=per_seq_gc_plot, width=26, height=2+2*nsamples, units="cm")
    
    per_seq_qual_plot = ggplot(data = per_seq_qual, aes(x=quality, y=norm_count)) +
        geom_col(fill="#114477") +
        scale_x_continuous(breaks=scales::pretty_breaks(n=5), name="quality score") +
        scale_y_continuous(breaks=scales::pretty_breaks(n=2), name="normalized counts") +
        facet_grid(sample~status, switch="y") +
        ggtitle("per sequence quality scores") +
        theme_standard +
        theme(axis.text.y = element_text(size=10, face="plain"))
    
    ggsave(per_seq_qual_out, plot=per_seq_qual_plot, width=26, height=2+1.5*nsamples, units="cm")
    
    dup_level_plot = ggplot(data = duplication_levels, aes(x=duplication_level, y=pct_of_total)) +
        geom_col(fill="#114477") +
        xlab("duplication level") +
        ylab("% of total reads") +
        facet_grid(sample~status, switch="y") +
        ggtitle("sequence duplication levels") +
        theme_standard +
        theme(axis.text.x = element_text(size=10, face="plain", angle=60, hjust=1),
              axis.text.y = element_text(size=10, face="plain"))
    
    ggsave(seq_dup_out, plot=dup_level_plot, width=26, height=2+1.5*nsamples, units="cm")
    
    #ermmm...no obvious way to make this one look nice, but then it doesn't really need to
    kmer_content_plot = ggplot(data = kmer_content, aes(x=max_position, y=log2(obs_over_exp_max), label=sequence)) +
        geom_point(shape=16, stroke=0, size=1, alpha=0.5) +
        geom_label_repel(size=2, label.size=unit(0.05, "pt"), label.padding=unit(0.1, "pt"), label.r=unit(0,"pt"), segment.size=0.1,
                         box.padding=unit(0.05,"pt"), segment.alpha=0.4) +
        xlab("position in read") +
        ylab(expression(bold(log[2]~ frac("observed", "expected")))) +
        ggtitle("k-mer content",
                subtitle = "top 20 overrepresented k-mers") +
        facet_grid(sample~status, switch="y", scales="free_y") +
        theme_standard + theme(plot.subtitle = element_text(size=12, face="plain"))
    
    ggsave(kmer_out, plot=kmer_content_plot, width=35, height=2+5*nsamples, units="cm", limitsize=FALSE)
}

main(seq_len_dist_in = snakemake@input[["seq_len_dist"]],
     per_tile_in = snakemake@input[["per_tile"]],
     per_base_qual_in = snakemake@input[["per_base_qual"]],
     per_base_seq_in = snakemake@input[["per_base_seq"]],
     per_base_n_in = snakemake@input[["per_base_n"]],
     per_seq_gc_in = snakemake@input[["per_seq_gc"]],
     per_seq_qual_in = snakemake@input[["per_seq_qual"]],
     adapter_content_in = snakemake@input[["adapter_content"]],
     seq_dup_in = snakemake@input[["seq_dup"]],
     kmer_in = snakemake@input[["kmer"]],
     seq_len_dist_out = snakemake@output[["seq_len_dist"]],
     per_tile_out = snakemake@output[["per_tile"]],
     per_base_qual_out = snakemake@output[["per_base_qual"]],
     per_base_seq_out = snakemake@output[["per_base_seq"]],
     per_seq_gc_out = snakemake@output[["per_seq_gc"]],
     per_seq_qual_out = snakemake@output[["per_seq_qual"]],
     adapter_content_out = snakemake@output[["adapter_content"]],
     seq_dup_out = snakemake@output[["seq_dup"]],
     kmer_out = snakemake@output[["kmer"]])
