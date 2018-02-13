library(tidyverse)
library(forcats)
library(gridExtra)
library(ggthemes)

main = function(intable, samplelist, controls, conditions, plotpath, statspath){
    df = read_table2(intable, col_names=c('sample', 'group', 'total','exp', 'si')) %>%
        filter(sample %in% samplelist) %>%
        mutate(sipct = si/total*100) %>%
        group_by(group) %>%
        mutate(outlier= ifelse(sipct>2.5*quantile(sipct,.75)-1.5*quantile(sipct,.25) |
                                   sipct< -2.5*quantile(sipct,.25)-1.5*quantile(sipct,.75),
                               TRUE, FALSE)) %>%
        ungroup() %>%
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))
    
    nsamples = nrow(df)
    ngroups = df$group %>% n_distinct()
    
    barplot = ggplot(data=df, aes(x=sample, fill=group, y=sipct)) +
        geom_col() +
        geom_text(aes(label=round(sipct, 1)), size=4,
                  position=position_stack(vjust=0.9)) +
        scale_fill_ptol(guide=FALSE) +
        ylab("% spike-in") +
        theme_bw() +
        theme(axis.text = element_text(size=12, color="black", face="bold"),
              axis.text.x = element_text(angle=30, hjust=0.9),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=12, color="black", face="bold",
                                          angle=0, vjust=0.5, hjust=1))
    
    boxplot = ggplot(data = df, aes(x=group, y=sipct, fill=group)) +
        geom_boxplot(outlier.shape=16, outlier.size=1.5, outlier.color="red", outlier.stroke=0) +
        geom_point(shape=16, size=1, stroke=0) +
        scale_fill_ptol(guide=FALSE) +
        ylab("% spike-in") +
        theme_bw() +
        theme(axis.text = element_text(size=12, color="black", face="bold"),
              axis.text.x = element_text(angle=30, hjust=0.9),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=12, color="black", face="bold"))
    
    outstats = df %>% add_count(group) %>% group_by(group) %>%
                mutate(median = median(sipct)) %>% ungroup() %>%
                filter(!outlier) %>% add_count(group) %>%
                group_by(group) %>%
                summarise(n = median(n), median = median(median),
                          n_no_outlier = median(nn),
                          mean_no_outlier = mean(sipct),
                          sd_no_outlier = sd(sipct)) 
    write_tsv(outstats, path = statspath, col_names=TRUE)
    
    ctrl = outstats %>% inner_join(tibble(c=controls), by=c("group"="c"))
    cond = outstats %>% inner_join(tibble(c=conditions), by=c("group"="c"))
    pct = bind_cols(ctrl, cond) %>%
            select(group, group1, mean_no_outlier, mean_no_outlier1) %>%
            rename(control=group, condition=group1,
                   pct.ctrl = mean_no_outlier, pct.cond = mean_no_outlier1) %>%
            mutate_at(c("pct.ctrl", "pct.cond"), funs(./100)) %>%
            mutate(relative.levels = pct.ctrl*(1-pct.cond)/(pct.cond*(1-pct.ctrl)))
    
    pctout = pct %>% select(control, condition, relative.levels) %>%
                mutate_at("relative.levels", funs(round(., digits=3)))
    pctdraw = tableGrob(pctout, rows=NULL, cols=c("control","condition","relative levels"),
                        ttheme_minimal(base_size=10))
    #set width
    wl = 1+1.4*nsamples
    wr = 1+1.6*ngroups
    th = 1+ngroups/4
    page = arrangeGrob(barplot, boxplot, pctdraw, layout_matrix=rbind(c(1,2),c(3,3)),
                       widths=unit(c(wl, wr), "cm"),
                       heights=unit(c(9,th),"cm"))
    
    ggsave(plotpath, page, width = wl+wr, height=9+th+.5, units = "cm")
}

df = main(intable = snakemake@input[[1]],
          samplelist = snakemake@params[["samplelist"]],
          controls = snakemake@params[["controls"]],
          conditions = snakemake@params[["conditions"]],
          plotpath = snakemake@output[["plot"]],
          statspath = snakemake@output[["stats"]])
