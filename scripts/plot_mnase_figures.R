library(tidyverse)
library(forcats)
library(magrittr)
library(viridis)
library(psych)
library(seriation)
library(ggthemes)
library(gtable)

#https://github.com/stas-g/findPeaks
find_peaks = function (x, m = 3){
    shape = diff(sign(diff(x, na.pad = FALSE)))
    pks = sapply(which(shape < 0), FUN = function(i){
        z = i - m + 1
        z = ifelse(z > 0, z, 1)
        w = i + m + 1
        w = ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
    pks = unlist(pks)
    pks
}

hmap_ybreaks = function(limits){
    if (max(limits)-min(limits) >= 2000){
        return(seq(min(limits)+500, max(limits)-500, 500))
    } else if (between(max(limits)-min(limits), 200, 1999)){
        return(seq(min(limits)+100, max(limits)-100, 100))
    } else {
        return(round((max(limits)-min(limits))/2))
    }
}


main = function(in_paths, samplelist, anno_paths, ptype, readtype, upstream, dnstream, scaled_length, pct_cutoff, spread_type,
                trim_pct, refptlabel, endlabel, cmap, sortmethod, cluster_scale, cluster_samples, cluster_five, cluster_three, k,
                heatmap_sample_out, heatmap_group_out, meta_sample_out, meta_sample_overlay_out, meta_group_out, meta_sampleclust_out, meta_groupclust_out,
                metahmap_sample_out, metahmap_group_out, anno_out, cluster_out){

    hmap = function(df, flimit){

        heatmap_base = ggplot(data = df) +
            geom_vline(xintercept = 0, size=1.5)

        if (ptype=="scaled"){
            heatmap_base = heatmap_base +
                geom_vline(xintercept = scaled_length/1000, size=1.5)
        }

        heatmap_base = heatmap_base +
            geom_raster(aes(x=position, y=new_index, fill=cpm), interpolate=FALSE) +
            scale_fill_viridis(option = cmap, na.value="#FFFFFF00",
                               name= if(readtype=="midpoint"){"MNase-seq dyad signal"} else {"MNase-seq protection"},
                               limits = c(NA, flimit), oob=scales::squish,
                               guide=guide_colorbar(title.position="top",
                                                    barwidth=20, barheight=1, title.hjust=0.5)) +
            scale_y_reverse(expand=c(0.005,5), breaks=hmap_ybreaks) +
            theme_minimal() +
            theme(text = element_text(size=16, face="bold", color="black"),
                  legend.position = "top",
                  legend.title = element_text(size=16, face="bold", color="black"),
                  legend.text = element_text(size=12, face="plain"),
                  legend.margin = margin(0,0,0,0),
                  legend.box.margin = margin(0,0,0,0),
                  strip.text.x = element_text(size=16, face="bold", color="black"),
                  axis.ticks.length = unit(0.125, "cm"),
                  axis.ticks = element_line(size=1.5),
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.text.x = element_text(size=16, face="bold", color="black", margin = unit(c(3,0,0,0),"pt")),
                  axis.title.x = element_text(size=12, face="plain"),
                  axis.title.y = element_blank(),
                  panel.grid.major.x = element_line(color="black", size=1.5),
                  panel.grid.minor.x = element_line(color="black"),
                  panel.grid.major.y = element_line(color="black"),
                  panel.grid.minor.y = element_blank(),
                  panel.spacing.x = unit(.8, "cm"))
        if (ptype=="absolute") {
            heatmap_base = heatmap_base +
                scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                                   labels= function(x){if_else(x==0, refptlabel,
                                                               if(upstream>500 | dnstream>500){as.character(x)}
                                                               else {as.character(x*1000)})},
                                   name=paste("distance from", refptlabel, if(upstream>500 | dnstream>500){"(kb)"}
                                              else {"(nt)"}),
                                   limits = c(-upstream/1000, dnstream/1000),
                                   expand=c(0,0.025))
        } else {
            heatmap_base = heatmap_base +
                scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels=c(refptlabel, "", endlabel),
                                   name="scaled distance",
                                   limits = c(-upstream/1000, (scaled_length+dnstream)/1000),
                                   expand=c(0,0.025))
        }
        return(heatmap_base)
    }

    meta = function(df, groupvar="sample", n){

        palette = if(n > 12){rep(ptol_pal()(12), ceiling(n/12))} else {ptol_pal()(n)}

        if (groupvar=="sample"){
            metagene = ggplot(data = df, aes(x=position, group=sample,
                                             color=group, fill=group))
        } else if (groupvar=="group"){
            metagene = ggplot(data = df, aes(x=position, group=group,
                                             color=group, fill=group))
        } else if (groupvar=="sampleclust"){
            metagene = ggplot(data = df %>% mutate(sampleclustid = paste(sample, cluster)),
                              aes(x=position, group=sampleclustid,
                                  color=factor(cluster), fill=factor(cluster)))
        } else if (groupvar=="groupclust"){
            metagene = ggplot(data = df %>% mutate(groupclustid = paste(group, cluster)),
                              aes(x=position, group=groupclustid,
                                  color=factor(cluster), fill=factor(cluster)))
        }

        metagene = metagene +
            geom_vline(xintercept = 0, size=1, color="grey65")

        if (ptype=="scaled"){
            metagene = metagene +
                geom_vline(xintercept = scaled_length/1000, size=1, color="grey65")
        }

        metagene = metagene +
            geom_ribbon(aes(ymin=low, ymax=high),
                        alpha=0.4, size=0) +
            geom_line(aes(y=mid)) +
            scale_y_continuous(limits = c(NA, NA), name="normalized counts") +
            scale_color_manual(values=palette,
                               guide=guide_legend(label.position="top", label.hjust=0.5)) +
            scale_fill_manual(values=palette) +
            ggtitle(if(readtype=="midpoint"){"MNase-seq dyad signal"} else {"MNase-seq protection"}) +
            theme_light() +
            theme(text = element_text(size=12, color="black", face="bold"),
                  axis.text = element_text(size=12, color="black"),
                  axis.text.y = element_text(size=10, face="plain"),
                  axis.title = element_text(size=10, face="plain"),
                  strip.placement="outside",
                  strip.background = element_blank(),
                  strip.text = element_text(size=12, color="black", face="bold"),
                  legend.text = element_text(size=12),
                  legend.title = element_blank(),
                  legend.position = "top",
                  legend.key.width = unit(3, "cm"),
                  plot.title = element_text(size=12),
                  plot.subtitle = element_text(size=10, face="plain"),
                  panel.spacing.x = unit(0.8, "cm"))
        if (ptype=="absolute"){
            metagene = metagene +
                scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                                   labels= function(x){if_else(x==0, refptlabel,
                                                               if(upstream>500 | dnstream>500){as.character(x)}
                                                               else {as.character(x*1000)})},
                                   name=paste("distance from", refptlabel, if(upstream>500 | dnstream>500){"(kb)"}
                                              else {"(nt)"}),
                                   limits = c(-upstream/1000, dnstream/1000),
                                   expand=c(0,0))
        } else {
            metagene = metagene +
                scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels=c(refptlabel, "", endlabel),
                                   name="scaled distance",
                                   limits = c(-upstream/1000, (scaled_length+dnstream)/1000),
                                   expand=c(0,0))
        }

        if (groupvar %in% c("sampleclust", "groupclust")){
            metagene = metagene +
                scale_color_colorblind(guide=guide_legend(label.position="top", label.hjust=0.5)) +
                scale_fill_colorblind()
        }
        return(metagene)
    }

    metahmap = function(meta_df, peak_df) {
        plot = ggplot() +
            geom_vline(xintercept = 0, size=1, color="grey65")
        if (ptype=="scaled"){
            plot = plot +
                geom_vline(xintercept = scaled_length, size=1, color="grey65") +
                scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels=c(refptlabel, "", endlabel),
                                   name="scaled distance",
                                   limits = c(-upstream/1000, (scaled_length+dnstream)/1000),
                                   expand=c(0,0))
        } else {
            plot = plot +
                scale_x_continuous(labels= function(x){if_else(x==0, refptlabel,
                                                               if(upstream>500 | dnstream>500){as.character(x)} else {as.character(x*1000)})},
                                   name=paste("distance from", refptlabel, if(upstream>500 | dnstream>500){"(kb)"} else {"(nt)"}),
                                   limits = c(-upstream/1000, dnstream/1000))
        }

        plot = plot +
            geom_vline(data = peak_df, aes(xintercept=position)) +
            geom_text(data = peak_df,
                      aes(x=position, y=-0.6,
                          label=if_else(position<=0, as.character(as.integer(position*1000)),
                                        paste0("+", as.integer(position*1000)))),
                      vjust=1, hjust=0, nudge_x=0.01, size=10/72*25.4) +
            geom_raster(data = meta_df, aes(x=position, y=0, fill=mid)) +
            scale_y_continuous(limits = c(-0.75, 0.5)) +
            scale_fill_viridis(option="inferno",
                               name = if(readtype=="midpoint"){"mean MNase dyad signal"} else {"mean MNase protection"},
                               guide=guide_colorbar(title.position="top", title.hjust=0.5, barwidth=15, barheight=1)) +
            theme_light() +
            theme(text = element_text(size=12, color="black", face="bold"),
                  strip.background = element_blank(),
                  strip.text = element_text(size=12, color="black"),
                  strip.text.y = element_text(angle=-180, hjust=1),
                  strip.placement="outside",
                  axis.ticks.y = element_blank(),
                  axis.text.x = element_text(size=12, color="black", face="bold"),
                  axis.text.y = element_blank(),
                  axis.title.x = element_text(size=10, face="plain"),
                  axis.title.y = element_blank(),
                  legend.text = element_text(size=10, face="plain"),
                  legend.position = "top",
                  legend.margin = margin(0,0,0,0),
                  legend.box.margin = margin(0,0,0,0),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.border = element_blank(),
                  panel.spacing.x = unit(0.5, "cm"))
        return(plot)
    }

    nest_right_facets = function(ggp, level=2, outer="replicate", inner="annotation"){
        og_grob = ggplotGrob(ggp)
        strip_loc = grep("strip-r", og_grob[["layout"]][["name"]])
        strip = gtable_filter(og_grob, "strip-r", trim=FALSE)
        strip_heights = gtable_filter(og_grob, "strip-r")[["heights"]]

        strip_top = min(strip[["layout"]][["t"]])
        strip_bot = max(strip[["layout"]][["b"]])
        strip_x = strip[["layout"]][["r"]][1]

        mat = matrix(vector("list", length=(length(strip)*2-1)*level), ncol=level)
        mat[] = list(zeroGrob())

        facet_grob = gtable_matrix("rightcol", grobs=mat,
                                   widths=unit(rep(1,level), "null"),
                                   heights=strip_heights)

        if(level==3){
            rep_grob_indices = seq(1, length(strip_loc), sum(k))
            for (rep_idx in 1:max_reps){
                #add replicate facet label
                facet_grob = facet_grob %>%
                        gtable_add_grob(grobs = og_grob$grobs[[strip_loc[rep_grob_indices[rep_idx]]]]$grobs[[level]],
                                        t = ((sum(k)*2))*(rep_idx-1)+1,
                                        b = ((sum(k)*2))*(rep_idx)-1,
                                        l = level, r = level)
                #for each annotation within each replicate
                for (anno_idx in 1:n_anno){
                    t = ((sum(k)*2))*(rep_idx-1)+1+sum(k[1:anno_idx])-k[1]+2*(anno_idx-1)
                    b = t + k[anno_idx]
                    facet_grob = facet_grob %>%
                            gtable_add_grob(grobs = og_grob$grobs[[strip_loc[rep_grob_indices[rep_idx]]+
                                                                       sum(k[1:anno_idx])-k[1]]]$grobs[[2]],
                                            t = t, b = b, l = 2, r = 2)
                }
            }
        } else if(level==2){
            if (outer=="annotation"){
                outer_grob_indices = 1+lag(k, default=0)
                n_outer = n_anno
            } else if (outer=="replicate"){
                outer_grob_indices = seq(1, length(strip_loc), sum(k))
                n_outer = max_reps
            }

            for (idx in 1:n_outer){
                if (outer=="annotation"){
                    t=((k[idx]*2))*(idx-1)+1
                    b=((k[idx]*2))*(idx)-1
                } else if (outer=="replicate"){
                    if (inner=="cluster"){
                        t=((k*2))*(idx-1)+1
                        b=((k*2))*(idx)-1
                    } else{
                        t = (n_anno*2)*(idx-1)+1
                        b = (n_anno*2)*(idx)-1
                    }
                }
                facet_grob %<>%
                    gtable_add_grob(grobs = og_grob$grobs[[strip_loc[outer_grob_indices[idx]]]]$grobs[[2]],
                                    t=t, b=b, l=2, r=2)
            }
        }
        new_grob = gtable_add_grob(og_grob, facet_grob, t=strip_top, r=strip_x, l=strip_x, b=strip_bot, name='rstrip')
        return(new_grob)
    }

    nest_top_facets = function(ggp, level=2, inner="cluster", intype="gg"){
        if (intype=="gg"){
            og_grob = ggplotGrob(ggp)
        } else if (intype=="gtable"){
            og_grob = ggp
        }

        strip_loc = grep("strip-t", og_grob[["layout"]][["name"]])
        strip = gtable_filter(og_grob, "strip-t", trim=FALSE)
        strip_widths = gtable_filter(og_grob, "strip-t")[["widths"]]

        strip_l = min(strip[["layout"]][["l"]])
        strip_r = max(strip[["layout"]][["r"]])
        strip_y = strip[["layout"]][["t"]][1]

        mat = matrix(vector("list", length=(length(strip)*2-1)*level), nrow=level)
        mat[] = list(zeroGrob())

        facet_grob = gtable_matrix("toprow", grobs=mat,
                                   heights=unit(rep(1,level), "null"),
                                   widths=strip_widths)
        if (inner=="cluster"){
            outer_grob_indices = 1+lag(k, default=0)
            n_outer = n_anno
        } else if (inner=="strand"){
            outer_grob_indices = seq(1, n_groups*2, 2)
            n_outer = n_groups
        }

        for (idx in 1:n_outer){
            if (inner=="cluster"){
                l=((k[idx]*2))*(idx-1)+1
                r=((k[idx]*2))*(idx)-1
            } else if (inner=="strand"){
                l=4*(idx-1)+1
                r=4*(idx)-1
            }
            facet_grob %<>%
                gtable_add_grob(grobs = og_grob$grobs[[strip_loc[outer_grob_indices[idx]]]]$grobs[[1]],
                                l=l, r=r, t=1, b=1)
        }
        new_grob = gtable_add_grob(og_grob, facet_grob, t=strip_y, r=strip_r, l=strip_l, b=strip_y, name='rstrip')
        return(new_grob)
    }

    import = function(path){
        read_tsv(path, col_names=c("group", "sample", "annotation", "index", "position","cpm")) %>%
            filter((sample %in% samplelist | sample %in% cluster_samples) & !is.na(cpm))
    }

    df = import(in_paths[1]) %>%
        group_by(annotation) %>%
        mutate(annotation_labeled = paste(n_distinct(index), annotation)) %>%
        ungroup() %>%
        mutate(annotation = annotation_labeled) %>%
        select(-annotation_labeled)%>%
        mutate_at(vars(group, sample, annotation), ~(fct_inorder(., ordered=TRUE)))

    #get replicate info for sample facetting
    repl_df = df %>%
        select(group, sample) %>%
        distinct() %>%
        group_by(group) %>%
        mutate(replicate=row_number()) %>% ungroup() %>% select(-group)
    max_reps = max(repl_df[["replicate"]])

    df %<>% left_join(repl_df, by="sample")

    n_anno = n_distinct(df[["annotation"]])

    #import annotation information
    annotations = df %>% distinct(annotation) %>% pull(annotation)
    bed = tibble()
    for (i in 1:n_anno){
        bed = read_tsv(anno_paths[i], col_names=c('chrom','start','end','name','score','strand')) %>%
            mutate(annotation=annotations[i]) %>%
            rowid_to_column(var="index") %>%
            bind_rows(bed, .)
    }

    n_samples = length(samplelist)
    n_groups = n_distinct(df[["group"]])

    #clustering
    if (sortmethod=="cluster"){
        reorder = tibble()

        #cluster for each annotation
        for (i in 1:length(annotations)){
            rr = df %>% filter(annotation==annotations[i] & sample %in% cluster_samples &
                                   between(position, cluster_five/1000, cluster_three/1000))
            if (cluster_scale){
                rr %<>%
                    group_by(group, sample, annotation, index, replicate) %>%
                    mutate(cpm = scales::rescale(cpm)) %>% ungroup()
            }
            rr %<>%
                select(-c(group, annotation, replicate)) %>%
                unite(cid, c(sample, position), sep="~") %>%
                spread(cid, cpm, fill=0) %>%
                select(-index)

            d = dist(rr, method="euclidean")
            l = kmeans(d, k[i])[["cluster"]]

            pdf(file=cluster_out[i], width=6, height=6)
            unsorted = dissplot(d, method=NA, newpage=TRUE,
                                main=paste0(annotations[i], "\nEuclidean distances, unsorted"),
                                options=list(silhouettes=FALSE, col=viridis(100, direction=-1)))
            if (k[i] > 1){
                seriated = dissplot(d, labels=l, method="OLO", newpage=TRUE,
                                    main=paste0(annotations[i], "\nEuclidean distances, ",
                                                k[i], "-means clustered,\nOLO inter- and intracluster sorting"),
                                    options=list(silhouettes=TRUE, col=viridis(100, direction=-1)))
                dev.off()

                sub_reorder = tibble(annotation = annotations[i],
                                     cluster = seriated[["labels"]],
                                     og_index = seriated[["order"]]) %>%
                    mutate(new_index = row_number())
            } else if (k[i]==1){
                seriated = seriate(d, method="OLO")
                sub_reorder = tibble(annotation = annotations[i],
                                     cluster = as.integer(1),
                                     og_index = get_order(seriated)) %>%
                    mutate(new_index = row_number())
            }
            reorder %<>% bind_rows(sub_reorder)

            sorted = sub_reorder %>%
                left_join(bed, by=c("annotation", "og_index"="index")) %>%
                select(-c(annotation, og_index, new_index))
            for (j in 1:k[i]){
                sorted %>% filter(cluster==j) %>%
                    select(-cluster) %>%
                    write_tsv(anno_out[sum(k[0:(i-1)])+j], col_names=FALSE)
            }
        }

        df %<>% left_join(reorder, by=c("annotation", "index"="og_index")) %>%
            group_by(annotation, cluster) %>%
            mutate(new_index = as.integer(new_index+1-min(new_index))) %>%
            ungroup() %>%
            arrange(annotation, cluster, new_index)
    } else if (sortmethod=="length"){
        sorted = bed %>% group_by(annotation) %>%
            arrange(end-start, .by_group=TRUE) %>%
            rowid_to_column(var= "new_index") %>%
            mutate(new_index = as.integer(new_index+1-min(new_index))) %>%
            ungroup()

        for (i in 1:n_anno){
            sorted %>% filter(annotation==annotations[i]) %>%
                select(-c(new_index, index, annotation)) %>%
                write_tsv(path=anno_out[i], col_names=FALSE)
        }

        df = sorted %>%
            select(index, new_index, annotation) %>%
            right_join(df, by=c("annotation", "index")) %>%
            mutate(cluster=as.integer(1))
    } else {
        df %<>% mutate(new_index = index,
                           cluster = as.integer(1))

        for (i in 1:n_anno){
            bed %>% filter(annotation==annotations[i]) %>%
                select(-c(index, annotation)) %>%
                write_tsv(path=anno_out[i], col_names=FALSE)
        }
    }

    df_sample = df %>%
        mutate(replicate = fct_inorder(paste("replicate", replicate), ordered=TRUE),
               cluster = fct_inorder(paste("cluster", cluster), ordered=TRUE))
    sample_cutoff= df_sample %>%
        filter(cpm > 0) %>%
        pull(cpm) %>%
        quantile(probs=pct_cutoff, na.rm=TRUE)

    df_group = df %>%
        group_by(group, annotation, position, cluster, new_index) %>%
        summarise(cpm = mean(cpm)) %>% ungroup() %>%
        mutate(cluster = fct_inorder(paste("cluster", cluster), ordered=TRUE))
    group_cutoff = df_group %>%
        filter(cpm > 0) %>%
        pull(cpm) %>%
        quantile(probs=pct_cutoff, na.rm=TRUE)

    # if the sortmethod isn't length, fill missing data with minimum signal
    # (only for heatmaps, don't want to influence metagene values)
    if (sortmethod != "length"){
        df_sample %<>%
            group_by(group, sample, annotation, replicate, cluster) %>%
            complete(new_index, position, fill=list(cpm=min(df_sample[["cpm"]]))) %>%
            ungroup()
        df_group %<>%
            group_by(group, annotation, cluster) %>%
            complete(new_index, position, fill=list(cpm=min(df_group[["cpm"]]))) %>%
            ungroup()
    }

    heatmap_sample = hmap(df_sample, sample_cutoff)
    heatmap_group = hmap(df_group, group_cutoff)
    if (n_anno==1 && max(k)==1){
        heatmap_sample = heatmap_sample +
            ylab(annotations[1]) +
            facet_grid(replicate ~ group, scales="free_y", space="free_y") +
            theme(axis.title.y = element_text(size=16, face="bold", color="black", angle=90),
                  strip.text.y = element_text(size=16, face="bold", color="black"))

        heatmap_group = heatmap_group +
            facet_grid(.~group) +
            ylab(annotations[1]) +
            theme(axis.title.y = element_text(size=16, face="bold", color="black", angle=90))
    } else if (n_anno==1 && max(k)>1){
        heatmap_sample = heatmap_sample +
            ylab(annotations[1]) +
            facet_grid(replicate + cluster ~ group, scales="free_y", space="free_y") +
            theme(axis.title.y = element_text(size=16, face="bold", color="black", angle=90),
                  strip.text.y = element_text(size=12, face="bold", color="black"),
                  strip.background = element_rect(fill="white", size=0))
        heatmap_sample %<>% nest_right_facets(level=2, outer="replicate", inner="cluster")

        heatmap_group = heatmap_group +
            ylab(annotations[1]) +
            facet_grid(cluster ~ group, scales="free_y", space="free_y") +
            theme(axis.title.y = element_text(size=16, face="bold", color="black", angle=90))
    } else if (n_anno>1 && max(k)==1){
        heatmap_sample = heatmap_sample +
            facet_grid(replicate + annotation ~ group, scales="free_y", space="free_y") +
            theme(strip.text.y = element_text(size=12, face="bold", color="black"),
                  strip.background = element_rect(fill="white", size=0))
        heatmap_sample %<>% nest_right_facets(level=2, outer="replicate")

        heatmap_group = heatmap_group +
            facet_grid(annotation ~ group, scales="free_y", space="free_y")
    } else if (n_anno>1 && max(k)>1){
        heatmap_sample = heatmap_sample +
            facet_grid(replicate + annotation + cluster ~ group,
                       scales="free_y", space="free_y") +
            theme(strip.text.y = element_text(size=12, face="bold", color="black"),
                  strip.background = element_rect(fill="white", size=0))
        heatmap_sample %<>% nest_right_facets(level=3)

        heatmap_group = heatmap_group +
            facet_grid(annotation + cluster ~ group, scales="free_y", space="free_y") +
            theme(strip.text.y = element_text(size=16, face="bold", color="black"),
                  strip.background = element_rect(fill="white", size=0))
        heatmap_group %<>% nest_right_facets(level=2, outer="annotation")
    }
    ggsave(heatmap_sample_out, plot=heatmap_sample, width=2+9*n_groups, height=10+10*max_reps, units="cm", limitsize=FALSE)
    ggsave(heatmap_group_out, plot=heatmap_group, width=2+9*n_groups, height=30, units="cm", limitsize=FALSE)

    metadf_sample = df %>%
        group_by(group, sample, annotation, position, cluster, replicate)

    if (spread_type=="conf_int"){
        metadf_sample %<>%
            summarise(mid = winsor.mean(cpm, trim=trim_pct),
                      sd = winsor.sd(cpm, trim=trim_pct)) %>%
            mutate(low = mid-sd,
                   high = mid+sd)

        #with SD correction for small sample sizes (Gurland and Tripathi 1971)
        metadf_group = metadf_sample %>%
            group_by(group, annotation, position, cluster) %>%
            summarise(sd = sd(mid),
                      n = n_distinct(replicate),
                      mid = mean(mid)) %>%
            mutate(sem = sqrt((n-1)/2)*gamma((n-1)/2)/gamma(n/2)*sd/sqrt(n))
    } else if (spread_type=="quantile") {

        metadf_sample %<>%
            summarise(mid = median(cpm),
                      low = quantile(cpm, probs=trim_pct),
                      high = quantile(cpm, probs=(1-trim_pct)))

        metadf_group = df %>%
            group_by(group, annotation, position, cluster) %>%
            summarise(mid = median(cpm),
                      low = quantile(cpm, probs=trim_pct),
                      high = quantile(cpm, probs=(1-trim_pct)))
    }

    metadf_sample %<>%
        ungroup() %>% arrange(replicate) %>%
        mutate(replicate = fct_inorder(paste("replicate", replicate), ordered=TRUE)) %>%
        arrange(cluster) %>%
        mutate(cluster = fct_inorder(paste("cluster", cluster), ordered=TRUE))

    peakdf_sample = metadf_sample %>%
        group_by(group, sample, annotation, cluster, replicate) %>%
        do(tibble(position = .$position,
                  fitted=loess(mid~position, data=., span=100/(upstream+dnstream))[["fitted"]])) %>%
        slice(find_peaks(fitted, m=10) %>% as.integer())

    metadf_group %<>%
        ungroup() %>%
        arrange(cluster) %>%
        mutate(cluster = fct_inorder(paste("cluster", cluster), ordered=TRUE))

    peakdf_group = metadf_group %>%
        group_by(group, annotation, cluster) %>%
        do(tibble(position = .$position,
                  fitted=loess(mid~position, data=., span=100/(upstream+dnstream))[["fitted"]])) %>%
        slice(find_peaks(fitted, m=10) %>% as.integer())

    meta_sample = meta(metadf_sample, n=n_groups)
    meta_group = meta(metadf_group, groupvar="group", n=n_groups)
    meta_sampleclust = meta(metadf_sample, groupvar="sampleclust", n=max(k))
    meta_groupclust = meta(metadf_group, groupvar="groupclust", n=max(k))
    metahmap_sample = metahmap(metadf_sample, peakdf_sample)
    metahmap_group = metahmap(metadf_group, peakdf_group)

    if (n_anno==1 && max(k)==1){
        meta_sample = meta_sample +
            scale_color_manual(values=rep("#4477AA", 100)) +
            scale_fill_manual(values=rep("#4477AA", 100)) +
            facet_grid(replicate~group) +
            ggtitle(if(readtype=="midpoint"){"MNase-seq dyad signal"} else {"MNase-seq protection"},
                    subtitle = annotations[1]) +
            theme(legend.position="none")

        palette = if(n_groups > 12){rep(ptol_pal()(12), ceiling(n_groups/12))} else {ptol_pal()(n_groups)}

        meta_sample_overlay = meta(metadf_sample, n=n_groups) +
            scale_color_manual(values=palette) +
            ggtitle(if(readtype=="midpoint"){"MNase-seq dyad signal"} else {"MNase-seq protection"},
                    subtitle = annotations[1]) +
            theme(legend.position="right",
                  legend.key.width=unit(0.8, "cm"))

        meta_group = meta_group +
            scale_color_manual(values=palette) +
            ggtitle(if(readtype=="midpoint"){"MNase-seq dyad signal"} else {"MNase-seq protection"},
                    subtitle = annotations[1]) +
            theme(legend.position="right",
                  legend.key.width=unit(0.8, "cm"))

        meta_sampleclust = meta_sampleclust +
            facet_grid(. ~ group) +
            ggtitle(if(readtype=="midpoint"){"MNase-seq dyad signal"} else {"MNase-seq protection"},
                    subtitle = annotations[1]) +
            theme(legend.position="right",
                  legend.key.width=unit(1, "cm"))

        meta_groupclust = meta_groupclust +
            facet_grid(. ~ group) +
            ggtitle(if(readtype=="midpoint"){"MNase-seq dyad signal"} else {"MNase-seq protection"},
                    subtitle = annotations[1]) +
            theme(legend.position="right",
                  legend.key.width=unit(1, "cm"))

        metahmap_sample = metahmap_sample +
            facet_grid(annotation+group~replicate, switch="y")

        metahmap_group = metahmap_group +
            facet_grid(annotation+group~., switch="y")

        ggsave(meta_sample_out, plot = meta_sample, width=3+7*n_groups, height=2+5*max_reps, units="cm", limitsize=FALSE)
        ggsave(meta_sample_overlay_out, plot = meta_sample_overlay, width=16, height=9, units="cm", limitsize=FALSE)
        ggsave(meta_group_out, plot = meta_group, width=16, height=9, units="cm", limitsize=FALSE)
        ggsave(meta_sampleclust_out, plot = meta_sampleclust, width=6+7*n_groups, height=10, units="cm", limitsize=FALSE)
        ggsave(meta_groupclust_out, plot = meta_groupclust, width=6+7*n_groups, height=10, units="cm", limitsize=FALSE)
        ggsave(metahmap_sample_out, plot = metahmap_sample, width=20+10*max_reps, height=4+4*n_anno, units="cm", limitsize=FALSE)
        ggsave(metahmap_group_out, plot = metahmap_group, width=30, height=4+4*n_anno, units="cm", limitsize=FALSE)
    } else if (n_anno>1 && max(k)==1){
        meta_sample = meta_sample + facet_grid(replicate ~ annotation)

        meta_sample_overlay = meta_sample + facet_grid(.~annotation)

        meta_group = meta_group + facet_grid(.~annotation)

        meta_sampleclust = meta_sampleclust + facet_grid(annotation ~ group)

        meta_groupclust = meta_groupclust + facet_grid(annotation ~ group)

        metahmap_sample = metahmap_sample + facet_grid(annotation+group~replicate, switch="y")

        metahmap_group = metahmap_group + facet_grid(annotation+group~., switch="y")
    } else if (max(k)>1){
        meta_sample = meta_sample +
            facet_grid(replicate ~ annotation + cluster) +
            theme(strip.background = element_rect(fill="white", size=0))
        meta_sample %<>% nest_top_facets(level=2)

        meta_sample_overlay = meta(metadf_sample, n=n_groups) +
            facet_grid(cluster ~ annotation)

        meta_group = meta_group + facet_grid(cluster ~ annotation)

        meta_sampleclust = meta_sampleclust +
            facet_grid(annotation ~ group) +
            theme(legend.key.width=unit(2, "cm"))

        meta_groupclust = meta_groupclust +
            facet_grid(annotation ~ group) +
            theme(legend.key.width=unit(2, "cm"))

        metahmap_sample = metahmap_sample + facet_grid(annotation+cluster+group~replicate, switch="y")

        metahmap_group = metahmap_group + facet_grid(annotation+cluster+group~., switch="y")
    }

    if (!(n_anno==1 && max(k)==1)){
        ggsave(meta_sample_out, plot = meta_sample, width=3+6*sum(k), height=2+5*max_reps, units="cm", limitsize=FALSE)
        ggsave(meta_sample_overlay_out, plot = meta_sample_overlay, width=3+7*n_anno, height=2+5*max(k), units="cm", limitsize=FALSE)
        ggsave(meta_group_out, plot = meta_group, width=3+7*n_anno, height=2+5*max(k), units="cm", limitsize=FALSE)
        ggsave(meta_sampleclust_out, plot = meta_sampleclust, width=3+7*n_groups, height=3+6*n_anno, units="cm", limitsize=FALSE)
        ggsave(meta_groupclust_out, plot = meta_groupclust, width=3+7*n_groups, height=3+6*n_anno, units="cm", limitsize=FALSE)
        ggsave(metahmap_sample_out, plot = metahmap_sample, width=12+12*max_reps, height=4+4*sum(k), units="cm", limitsize=FALSE)
        ggsave(metahmap_group_out, plot = metahmap_group, width=30, height=4+4*sum(k), units="cm", limitsize=FALSE)
    }
}

main(in_paths = snakemake@input[["matrix"]],
     samplelist = snakemake@params[["samplelist"]],
     anno_paths = snakemake@input[["annotations"]],
     ptype = snakemake@params[["plottype"]],
     readtype = snakemake@wildcards[["readtype"]],
     upstream = snakemake@params[["upstream"]],
     dnstream = snakemake@params[["dnstream"]],
     scaled_length = snakemake@params[["scaled_length"]],
     pct_cutoff = snakemake@params[["pct_cutoff"]],
     spread_type = snakemake@params[["spread_type"]],
     trim_pct = snakemake@params[["trim_pct"]],
     refptlabel = snakemake@params[["refpointlabel"]],
     endlabel = snakemake@params[["endlabel"]],
     cmap = snakemake@params[["cmap"]],
     sortmethod = snakemake@params[["sortmethod"]],
     cluster_scale = snakemake@params[["cluster_scale"]],
     cluster_samples = snakemake@params[["cluster_samples"]],
     cluster_five = snakemake@params[["cluster_five"]],
     cluster_three = snakemake@params[["cluster_three"]],
     k = snakemake@params[["k"]],
     heatmap_sample_out = snakemake@output[["heatmap_sample"]],
     heatmap_group_out = snakemake@output[["heatmap_group"]],
     meta_sample_out = snakemake@output[["meta_sample"]],
     meta_sample_overlay_out = snakemake@output[["meta_sample_overlay"]],
     meta_group_out = snakemake@output[["meta_group"]],
     meta_sampleclust_out = snakemake@output[["meta_sample_clust"]],
     meta_groupclust_out = snakemake@output[["meta_group_clust"]],
     metahmap_sample_out = snakemake@output[["metahmap_sample"]],
     metahmap_group_out = snakemake@output[["metahmap_group"]],
     anno_out = snakemake@params[["annotations_out"]],
     cluster_out = snakemake@params[["cluster_out"]])

