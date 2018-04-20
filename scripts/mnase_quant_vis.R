library(tidyverse)
library(forcats)
library(viridis)
library(psych)
library(ggthemes)
library(pals)

main = function(individual_path, integrated_path, anno_paths, anno_labels, sortmethod, refpt, refptlabel, condition, control, upstream,
                max_length, trim_pct, binsize, occupancy_cutoffs, fuzziness_cutoffs, occupancy_lfc_limit, fuzziness_lfc_limit, displacement_limit,
                indiv_occ_hmap_out, indiv_fuzz_hmap_out, indiv_occ_meta_out, indiv_fuzz_meta_out,
                integrated_occ_summit_hmap_out, integrated_occ_point_hmap_out, integrated_fuzz_hmap_out, integrated_displacement_hmap_out, integrated_displacement_segment_hmap_out,
                integrated_occ_summit_meta_out, integrated_occ_point_meta_out, integrated_fuzz_meta_out, integrated_displacement_meta_out){
    #####
    # import data for individual conditions
    n_anno = length(anno_paths)
    annotation = tibble()
    for (i in 1:n_anno){
        annotation = read_tsv(anno_paths[i], col_names=c("chrom", "start", "end", "feat_name", "score", "feat_strand")) %>% 
            select(-score) %>%
            mutate(annotation=anno_labels[i]) %>% 
            bind_rows(annotation, .)
    }
    
    individual = read_tsv(individual_path) %>% 
        left_join(annotation, by=c("chrom", "feat_name", "feat_strand", "annotation")) %>% 
        mutate(feat_start=start, feat_end=end) %>% 
        select(-c(start, end))
    if(refpt=="center"){
        individual = individual %>%
            mutate(start = as.integer((start+end)/2),
                   end = as.integer(start+1))
    }
    individual = individual %>%
        mutate_at(vars(nuc_start, nuc_end, nuc_summit),
                  funs(if_else(feat_strand=="+", .-feat_start, feat_end-.))) %>% 
        group_by(annotation) %>% 
        mutate(anno_labeled = paste(n_distinct(feat_name), annotation)) %>% 
        ungroup() %>% mutate(annotation=anno_labeled) %>% select(-anno_labeled) %>%  
        mutate_at(vars(group, annotation),
                  funs(fct_inorder(., ordered=TRUE)))
    if (sortmethod=="length"){
        individual = individual %>% 
            group_by(annotation, feat_name) %>%
            arrange(annotation, feat_end-feat_start) %>% 
            ungroup() %>% 
            mutate(feat_name = fct_inorder(feat_name, ordered=TRUE))
    } else {
        individual = individual %>% 
            mutate(feat_name=fct_inorder(feat_name, ordered=TRUE))
    }
    
    individual_meta = individual %>% 
        mutate(position = cut(nuc_summit, breaks=seq(min(nuc_summit, na.rm=TRUE), max(nuc_summit, na.rm=TRUE), by=binsize),
                              labels=seq(min(nuc_summit, na.rm=TRUE)+binsize/2, max(nuc_summit, na.rm=TRUE)-binsize/2, by=binsize)) %>%
                   as.character() %>% as.integer()) %>% 
        filter(between(position, -upstream, max_length)) %>% 
        group_by(group, position, annotation) %>% 
        summarise(occupancy_mean = winsor.mean(occupancy, trim=trim_pct),
                  occupancy_sem = winsor.sd(occupancy, trim=trim_pct)/sqrt(n()),
                  fuzziness_mean = winsor.mean(fuzziness, trim=trim_pct),
                  fuzziness_sem = winsor.sd(fuzziness, trim=trim_pct)/sqrt(n())) %>% 
        drop_na()
    
    #####
    # import data for DANPOS2 integrated analysis
    integrated = read_tsv(integrated_path,
                          col_types = "ciicdcciiiiiiidddddddddddddddic") %>% 
        left_join(annotation, by=c("feat_chrom"="chrom", "feat_name", "feat_strand", "annotation")) %>% 
        mutate(feat_start=start, feat_end=end) %>% 
        select(-c(start, end, nuc_chrom, overlap)) %>% 
        mutate_at(vars(nuc_start, nuc_end, nuc_center, ctrl_summit_loc, cond_summit_loc, diff_summit_loc),
                  funs(if_else(feat_strand=="+", .-feat_start, feat_end-.))) %>% 
        group_by(annotation) %>% 
        mutate(anno_labeled = paste(n_distinct(feat_name), annotation)) %>% 
        ungroup() %>% mutate(annotation=anno_labeled) %>% select(-anno_labeled) %>%  
        mutate(cond_ctrl_dist = cond_summit_loc-ctrl_summit_loc,
               annotation = fct_inorder(annotation, ordered=TRUE)) %>% 
        mutate(direction=factor(as.integer(sign(cond_ctrl_dist)),
                                levels=c(-1, 0, 1),
                                labels=c("-", "no change", "+")))
    if (sortmethod=="length") {
        integrated = integrated %>% 
            group_by(feat_name) %>% 
            arrange(annotation, feat_end-feat_start) %>% 
            ungroup() %>% 
            mutate(feat_name = fct_inorder(feat_name, ordered=TRUE))
    } else {
        integrated = integrated %>% 
            mutate(feat_name = fct_inorder(feat_name, ordered=TRUE))
    }
    
    integrated_meta = integrated %>% 
        mutate(position = cut(nuc_center, breaks=seq(min(nuc_center, na.rm=TRUE), max(nuc_center, na.rm=TRUE), by=binsize),
                              labels=seq(min(nuc_center, na.rm=TRUE)+binsize/2, max(nuc_center, na.rm=TRUE)-binsize/2, by=binsize)) %>% 
                   as.character() %>% as.integer()) %>% 
        filter(between(position, -upstream, max_length)) %>% 
        group_by(annotation, position) %>% 
        summarise(summit_lfc_mean = winsor.mean(summit_lfc, trim=trim_pct),
                  summit_lfc_sem = winsor.sd(summit_lfc, trim=trim_pct)/sqrt(n()),
                  point_lfc_mean = winsor.mean(point_lfc, trim=trim_pct),
                  point_lfc_sem = winsor.sd(point_lfc, trim=trim_pct)/sqrt(n()),
                  fuzziness_lfc_mean = winsor.mean(fuzziness_lfc, trim=trim_pct),
                  fuzziness_lfc_sem = winsor.sd(fuzziness_lfc, trim=trim_pct)/sqrt(n()),
                  cond_ctrl_dist_mean = winsor.mean(cond_ctrl_dist, trim=trim_pct),
                  cond_ctrl_dist_sem = winsor.sd(cond_ctrl_dist, trim=trim_pct)/sqrt(n())) %>% 
        drop_na()
    
    #####
    # individual plots
    theme_indiv_heatmap = theme_dark() +
        theme(text = element_text(size=16, color="black", face="bold"),
              strip.background = element_blank(),
              strip.text = element_text(size=16, color="black", face="bold"),
              strip.text.y = if(n_anno==1){element_text(angle=-90)} else{element_text(angle=-180, hjust=1, vjust=0.5)},
              axis.text.x = element_text(size=16, color="black", face="bold"),
              axis.text.y = element_blank(),
              axis.title.x = element_text(size=12, face="plain"),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.background = element_rect(fill="black"),
              legend.position = "top",
              legend.text = element_text(size=12, face="plain"),
              legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(0,0,0,0),
              panel.spacing.x = unit(1, "cm"),
              plot.margin = margin(5.5, 16, 5.5, 5.5, "pt"))
    
    xscale = scale_x_continuous(limits = c(-upstream, max_length),
                                expand = c(0,0),
                                name=paste("distance from", refptlabel, ifelse(max_length>500,"(kb)","(nt)")),
                                labels = function(x)ifelse(x==0, refptlabel, if(max_length>500){x/1000}else{x}))
    
    indiv_occ_hmap = ggplot(data = individual,
           aes(x=nuc_summit, y=as.integer(feat_name), fill=occupancy)) +
        geom_tile(width=100, height=1, size=0, na.rm=TRUE) +
        xscale +
        scale_y_reverse(expand=c(0,0)) +
        scale_fill_viridis(option="inferno", oob=scales::squish, 
                           limits= individual %>% pull(occupancy) %>% quantile(probs=occupancy_cutoffs),
                           labels = scales::scientific,
                           name = "nucleosome occupancy",
                           guide=guide_colorbar(title.position="top", title.hjust=0.5,
                                                barwidth=16, barheight=1)) +
        facet_grid(annotation~group, switch="y", scales="free_y", space="free_y") +
        theme_indiv_heatmap
    
    indiv_fuzz_hmap = ggplot(data = individual,
           aes(x=nuc_summit, y=as.integer(feat_name), fill=fuzziness)) +
        geom_tile(width=100, height=1, size=0, na.rm=TRUE) +
        xscale +
        scale_y_reverse(expand=c(0,0)) +
        scale_fill_viridis(option="inferno", oob=scales::squish, 
                           limits= individual %>% pull(fuzziness) %>% quantile(probs=fuzziness_cutoffs),
                           name = "nucleosome fuzziness (bp)",
                           guide=guide_colorbar(title.position="top", title.hjust=0.5,
                                                barwidth=16, barheight=1)) +
        facet_grid(annotation~group, switch="y", scales="free_y", space="free_y") +
        theme_indiv_heatmap
    
    theme_indiv_meta = theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              axis.text.y = element_text(size=10, face="plain"),
              axis.title.x = element_text(size=10, face="plain"),
              legend.title = element_blank(),
              legend.text = element_text(size=12),
              legend.position= ifelse(n_anno==1, "right", "top"),
              strip.placement = "outside",
              strip.background = element_blank(),
              strip.text = element_text(size=12, color="black", face="bold"),
              strip.text.y = if(n_anno==1){element_blank()}else{element_text(angle=0, hjust=0, vjust=0.5)},
              plot.title = element_text(size=12),
              plot.subtitle = element_text(size=12, face="plain"))
    
    indiv_meta_func = function(ggp){
        ggp +
            geom_vline(xintercept = 0, size=1, color="grey65") +
            geom_ribbon(alpha=0.4, size=0) +
            geom_line() +
            xscale +
            scale_fill_ptol(guide=guide_legend(keywidth=unit(ifelse(n_anno==1,0.75,3), "cm"),
                                               label.position=ifelse(n_anno==1,"right","top"),
                                               label.hjust=ifelse(n_anno==1,0,0.5))) +
            scale_color_ptol() +
            facet_grid(annotation~., switch=if(n_anno==1){"y"}else{NULL}) +
            theme_indiv_meta
    }
    
    indiv_occ_meta = ggplot(data = individual_meta,
           aes(x=position, y=occupancy_mean, ymin=occupancy_mean-1.96*occupancy_sem,
               ymax=occupancy_mean+1.96*occupancy_sem, color=group, fill=group)) %>% indiv_meta_func() +
        scale_y_continuous(name="occupancy",
                           labels = scales::scientific) +
        ggtitle(label="mean nucleosome occupancy",
                subtitle = if(n_anno==1){levels(individual_meta[["annotation"]])[1]}else{NULL})
    
    indiv_fuzz_meta = ggplot(data = individual_meta,
           aes(x=position, y=fuzziness_mean, ymin=fuzziness_mean-1.96*fuzziness_sem,
               ymax=fuzziness_mean+1.96*fuzziness_sem, color=group, fill=group)) %>% indiv_meta_func() +
        scale_y_continuous(name=expression(bold(fuzziness %==% std. ~ dev. ~ of ~ positions ~ (bp)))) +
        ggtitle(label="mean nucleosome fuzziness",
                subtitle = if(n_anno==1){levels(individual_meta[["annotation"]])[1]}else{NULL}) +
        theme(axis.text.y = element_text(size=12, face="bold"))
    
    #####
    # integrated plots
    theme_integrated_heatmap = theme_light() +
        theme(text = element_text(size=16, color="black", face="bold"),
              axis.text = element_text(size=16, color="black"),
              axis.text.y = element_blank(),
              axis.title.x = element_text(size=12, face="plain"),
              axis.title.y = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size=16, color="black"),
              strip.text.y = if(n_anno==1){element_text(angle=-90)} else{element_text(angle=0, hjust=0)},
              strip.placement="outside",
              panel.border = element_blank(),
              panel.grid.major.x = element_line(color="grey65", size=1),
              panel.grid.minor.x = element_line(color="grey65"),
              axis.ticks.x = element_line(color="grey65", size=1),
              legend.position="top",
              legend.title = element_text(size=12, face="plain"),
              legend.text = element_text(size=12, face="plain"),
              legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(0,0,0,0),
              plot.title = element_text(size=16))
    
    integrated_hmap_func = function(ggp, lowerlimit, upperlimit){
        ggp +
            geom_tile(width=100, height=1, na.rm=TRUE) +
            scale_y_reverse(expand=c(0,0)) +
            xscale +
            scale_fill_gradientn(colors=coolwarm(100),
                                 limits=c(lowerlimit,upperlimit), oob=scales::squish,
                                 name=bquote(log[2] ~ frac(.(condition), .(control))),
                                 guide=guide_colorbar(title.position="left", title.vjust=1,
                                                      barwidth=16, barheight=1)) +
            facet_grid(annotation~., switch=if(n_anno==1){"y"}else{NULL}, scales="free_y", space="free_y") +
            theme_integrated_heatmap
    }
    
    integrated_occ_summit_hmap = ggplot(data = integrated,
           aes(x=nuc_center, y=as.integer(feat_name), fill=summit_lfc)) %>%
        integrated_hmap_func(lowerlimit=-occupancy_lfc_limit, upperlimit=occupancy_lfc_limit) +
        ggtitle("nucleosome occupancy changes (summit)") 
    
    integrated_occ_point_hmap = ggplot(data = integrated,
           aes(x=nuc_center, y=as.integer(feat_name), fill=point_lfc)) %>% 
        integrated_hmap_func(lowerlimit=-occupancy_lfc_limit, upperlimit=occupancy_lfc_limit) +
        ggtitle("nucleosome occupancy changes (point)")
    
    integrated_fuzz_hmap = ggplot(data = integrated, aes(x=nuc_center, y=as.integer(feat_name), fill=fuzziness_lfc)) %>% 
        integrated_hmap_func(lowerlimit=-fuzziness_lfc_limit, upperlimit=fuzziness_lfc_limit) +
        ggtitle("nucleosome fuzziness changes")
    
    integrated_displacement_hmap = ggplot(data = integrated, aes(x=nuc_center, y=as.integer(feat_name), fill=cond_ctrl_dist)) +
        geom_tile(width=100, height=1, na.rm=TRUE) +
        scale_y_reverse(expand=c(0,0)) +
        xscale +
        scale_fill_gradientn(colors=coolwarm(100), limits=c(-displacement_limit,displacement_limit), oob=scales::squish,
                             breaks = scales::pretty_breaks(5),
                             labels = function(x){if_else(x<=0, as.character(x), paste0("+",x))},
                             name="nucleosome displacement (bp)",
                             guide=guide_colorbar(title.position="top", title.hjust=0.5,
                                                  barwidth=16, barheight=1)) +
        facet_grid(annotation~., switch=if(n_anno==1){"y"}else{NULL}, scales="free_y", space="free_y") +
        theme_integrated_heatmap +
        theme(legend.title = element_text(size=16, face="bold"),
              legend.text = element_text(size=12, face="bold"))
    
    integrated_displacement_segment_hmap = ggplot(data = integrated, aes(x=ctrl_summit_loc, y=as.integer(feat_name),
                                  xend=cond_summit_loc, yend=as.integer(feat_name),
                                  color=direction)) +
        geom_segment(na.rm=TRUE) +
        scale_y_reverse(expand=c(0,0)) +
        xscale +
        scale_color_ptol(name="direction of nucleosome displacement:",
                         guide=guide_legend(title.position="top", label.position="top",
                                            label.hjust=0.5, keywidth=unit(2, "cm"))) +
        facet_grid(annotation~., switch=if(n_anno==1){"y"}else{NULL}, scales="free_y", space="free_y") +
        theme_integrated_heatmap +
        theme(legend.title = element_text(size=16, face="bold"),
              legend.text = element_text(size=12, face="bold"))
    
    theme_integrated_meta = theme_indiv_meta +
        theme(axis.text.y = element_text(size=12, face="bold"),
              axis.title.y = element_text(angle=0, vjust=0.5),
              legend.position = ifelse(n_anno==1, "none", "right"))
    
    integrated_meta_func = function(ggp){
        ggp +
            geom_hline(yintercept = 0, size=1, color="gray65") +
            geom_vline(xintercept = 0, size=1, color="grey65") +
            geom_ribbon(size=0, alpha=ifelse(n_anno<=3,0.4,0.1)) +
            geom_line(alpha=ifelse(n_anno<=3, 1, 0.5)) +
            xscale +
            scale_color_ptol() +
            scale_fill_ptol() +
            ylab(bquote(bold(log[2] ~ frac(.(condition), .(control))))) +
            theme_integrated_meta
    }
    
    integrated_occ_summit_meta = ggplot(data = integrated_meta,
                                        aes(x=position, y=summit_lfc_mean,
                                            ymin=summit_lfc_mean-1.96*summit_lfc_sem,
                                            ymax=summit_lfc_mean+1.96*summit_lfc_sem,
                                            color=annotation, fill=annotation)) %>%
        integrated_meta_func()+
        ggtitle(label="mean nucleosome occupancy changes (summit)",
                subtitle = if(n_anno==1){levels(individual_meta[["annotation"]])[1]}else{NULL})
    
    integrated_occ_point_meta = ggplot(data = integrated_meta,
                                       aes(x=position, y=point_lfc_mean,
                                           ymin=point_lfc_mean-1.96*point_lfc_sem,
                                           ymax=point_lfc_mean+1.96*point_lfc_sem,
                                           color=annotation, fill=annotation)) %>% 
        integrated_meta_func()+
        ggtitle(label="mean nucleosome occupancy changes (point)",
                subtitle = if(n_anno==1){levels(individual_meta[["annotation"]])[1]}else{NULL})
    
    integrated_fuzz_meta = ggplot(data = integrated_meta,
                                  aes(x=position, y=fuzziness_lfc_mean,
                                      ymin=fuzziness_lfc_mean-1.96*fuzziness_lfc_sem,
                                      ymax=fuzziness_lfc_mean+1.96*fuzziness_lfc_sem,
                                      color=annotation, fill=annotation)) %>% 
        integrated_meta_func()+
        ggtitle(label="mean nucleosome fuzziness changes",
                subtitle = if(n_anno==1){levels(individual_meta[["annotation"]])[1]}else{NULL})
    
    integrated_displacement_meta = ggplot(data = integrated_meta,
                                          aes(x=position, y=cond_ctrl_dist_mean,
                                              ymin=cond_ctrl_dist_mean-1.96*cond_ctrl_dist_sem,
                                              ymax=cond_ctrl_dist_mean+1.96*cond_ctrl_dist_sem,
                                              color=annotation, fill=annotation)) +
        geom_hline(yintercept = 0, size=1, color="gray65") +
        geom_vline(xintercept = 0, size=1, color="grey65") +
        geom_ribbon(size=0, alpha=ifelse(n_anno<=3,0.4,0.1)) +
        geom_line(alpha=ifelse(n_anno<=3, 1, 0.5)) +
        xscale +
        scale_color_ptol() +
        scale_fill_ptol() +
        scale_y_continuous(breaks = scales::pretty_breaks(n=8),
                           labels = function(x){if_else(x<=0, as.character(x), paste0("+",x))},
                           name="displacement (bp)") +
        ggtitle(label="mean nucleosome displacement",
                subtitle = if(n_anno==1){levels(individual_meta[["annotation"]])[1]}else{NULL}) +
        theme_integrated_meta
    
    ggsave(indiv_occ_hmap_out, plot=indiv_occ_hmap, width=ifelse(n_anno==1,20,25), height=30, units="cm")
    ggsave(indiv_fuzz_hmap_out, plot=indiv_fuzz_hmap, width=ifelse(n_anno==1,20,25), height=30, units="cm")
    ggsave(indiv_occ_meta_out, plot=indiv_occ_meta, width=ifelse(n_anno==1,16,22), height=2+8*n_anno, units="cm")
    ggsave(indiv_fuzz_meta_out, plot=indiv_fuzz_meta, width=ifelse(n_anno==1,16,22), height=2+8*n_anno, units="cm")
    ggsave(integrated_occ_summit_hmap_out, plot=integrated_occ_summit_hmap, width=ifelse(n_anno==1,16,22), height=30, units="cm")
    ggsave(integrated_occ_point_hmap_out, plot=integrated_occ_point_hmap, width=ifelse(n_anno==1,16,22), height=30, units="cm")
    ggsave(integrated_fuzz_hmap_out, plot=integrated_fuzz_hmap, width=ifelse(n_anno==1,16,22), height=30, units="cm")
    ggsave(integrated_displacement_hmap_out, plot=integrated_displacement_hmap, width=ifelse(n_anno==1,16,22), height=30, units="cm")
    ggsave(integrated_displacement_segment_hmap_out, plot=integrated_displacement_segment_hmap, width=ifelse(n_anno==1,16,22), height=30, units="cm")
    ggsave(integrated_occ_summit_meta_out, plot=integrated_occ_summit_meta, width=ifelse(n_anno==1,16,22), height=10, units="cm")
    ggsave(integrated_occ_point_meta_out, plot=integrated_occ_point_meta, width=ifelse(n_anno==1,16,22), height=10, units="cm")
    ggsave(integrated_fuzz_meta_out, plot=integrated_fuzz_meta, width=ifelse(n_anno==1,16,22), height=10, units="cm")
    ggsave(integrated_displacement_meta_out, plot=integrated_displacement_meta, width=ifelse(n_anno==1,16,22), height=10, units="cm")

}

main(sortmethod= snakemake@params[["sortmethod"]],
     refpt = snakemake@params[["refpoint"]],
     refptlabel = snakemake@params[["refptlabel"]],
     condition= snakemake@wildcards[["condition"]],
     control= snakemake@wildcards[["control"]],
     individual_path = snakemake@input[["individual"]],
     integrated_path = snakemake@input[["integrated"]],
     anno_paths = snakemake@input[["annotations"]],
     anno_labels = snakemake@params[["anno_labels"]],
     binsize= snakemake@params[["binsize"]],
     occupancy_cutoffs = snakemake@params[["occupancy_cutoffs"]],
     fuzziness_cutoffs = snakemake@params[["fuzziness_cutoffs"]],
     fuzziness_lfc_limit = snakemake@params[["fuzziness_lfc_limit"]],
     displacement_limit = snakemake@params[["displacement_limit"]],
     occupancy_lfc_limit = snakemake@params[["occupancy_lfc_limit"]],
     upstream = snakemake@params[["upstream"]],
     max_length = snakemake@params[["max_length"]],
     trim_pct = snakemake@params[["trim_pct"]],
     indiv_occ_hmap_out = snakemake@output[["indiv_occ_hmap"]],
     indiv_fuzz_hmap_out = snakemake@output[["indiv_fuzz_hmap"]],
     indiv_occ_meta_out = snakemake@output[["indiv_occ_meta"]],
     indiv_fuzz_meta_out = snakemake@output[["indiv_fuzz_meta"]],
     integrated_occ_summit_hmap_out = snakemake@output[["integrated_occ_summit_hmap"]],
     integrated_occ_point_hmap_out = snakemake@output[["integrated_occ_point_hmap"]],
     integrated_fuzz_hmap_out = snakemake@output[["integrated_fuzz_hmap"]],
     integrated_displacement_hmap_out = snakemake@output[["integrated_displacement_hmap"]],
     integrated_displacement_segment_hmap_out = snakemake@output[["integrated_displacement_segment_hmap"]],
     integrated_occ_summit_meta_out = snakemake@output[["integrated_occ_summit_meta"]],
     integrated_occ_point_meta_out = snakemake@output[["integrated_occ_point_meta"]],
     integrated_fuzz_meta_out = snakemake@output[["integrated_fuzz_meta"]],
     integrated_displacement_meta_out = snakemake@output[["integrated_displacement_meta"]])
