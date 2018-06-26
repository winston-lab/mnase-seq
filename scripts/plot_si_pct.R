library(tidyverse)
library(forcats)
library(gridExtra)
library(ggthemes)

main = function(in_path, sample_list, controls, conditions, plot_out, stats_out){
    df = read_tsv(in_path) %>%
        filter(sample %in% sample_list) %>%
        mutate(si_pct = spikein_counts/total_counts*100) %>%
        group_by(group) %>%
        mutate(outlier= ifelse(si_pct > 2.5*quantile(si_pct,.75) - 1.5*quantile(si_pct,.25) |
                                   si_pct< -2.5*quantile(si_pct,.25) - 1.5*quantile(si_pct,.75),
                               TRUE, FALSE)) %>%
        ungroup() %>%
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))

    n_samples = nrow(df)
    n_groups = df %>% pull(group) %>% n_distinct()

    barplot = ggplot(data=df, aes(x=sample, fill=group, y=si_pct)) +
        geom_col() +
        geom_text(aes(label=round(si_pct, 1)), size=12/75*25.4,
                  position=position_stack(vjust=0.9)) +
        scale_fill_ptol(guide=FALSE) +
        ylab("% spike-in") +
        theme_light() +
        theme(axis.text = element_text(size=12, color="black", face="bold"),
              axis.text.x = element_text(angle=30, hjust=0.9),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=12, color="black", face="bold",
                                          angle=0, vjust=0.5, hjust=1))

    boxplot = ggplot(data = df, aes(x=group, y=si_pct, fill=group)) +
        geom_boxplot(outlier.shape=16, outlier.size=1.5, outlier.color="red", outlier.stroke=0) +
        geom_point(shape=16, size=1, stroke=0) +
        scale_fill_ptol(guide=FALSE) +
        ylab("% spike-in") +
        theme_light() +
        theme(axis.text = element_text(size=12, color="black", face="bold"),
              axis.text.x = element_text(angle=30, hjust=0.9),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=12, color="black", face="bold"))

    si_stats = df %>%
        add_count(group) %>%
        group_by(group) %>%
        mutate(median = median(si_pct)) %>%
        ungroup() %>%
        filter(!outlier) %>%
        add_count(group) %>%
        group_by(group) %>%
        summarise(n = median(n), median = median(median),
                  n_no_outlier = median(nn),
                  mean_no_outlier = mean(si_pct),
                  sd_no_outlier = sd(si_pct)) %>%
        write_tsv(path = stats_out, col_names=TRUE)

    rel_levels = tibble(condition=conditions, control=controls) %>%
        left_join(si_stats %>% select(group, mean_no_outlier),
                  by=c("condition"="group")) %>%
        rename(condition_pct=mean_no_outlier) %>%
        left_join(si_stats %>% select(group, mean_no_outlier),
                  by=c("control"="group")) %>%
        rename(control_pct=mean_no_outlier) %>%
        mutate_at(vars(condition_pct, control_pct),
                  funs(./100)) %>%
        mutate(relative_levels = control_pct*(1-condition_pct)/(condition_pct*(1-control_pct)))

    pct_table = rel_levels %>%
        select(condition, control, relative_levels) %>%
        mutate_at("relative_levels", funs(round(., digits=3)))
    pct_draw = tableGrob(pct_table, rows=NULL, cols=c("condition","control","relative levels"),
                        ttheme_minimal(base_size=10))
    #set width
    wl = 1+1.6*n_samples
    wr = 1+1.8*n_groups
    th = 1+length(conditions)/2
    page = arrangeGrob(barplot, boxplot, pct_draw, layout_matrix=rbind(c(1,2),c(3,3)),
                       widths=unit(c(wl, wr), "cm"),
                       heights=unit(c(9,th),"cm"))

    ggsave(plot_out, page, width = wl+wr, height=9+th+.5, units = "cm")
}

main(in_path = snakemake@input[[1]],
     sample_list = snakemake@params[["samplelist"]],
     controls = snakemake@params[["controls"]],
     conditions = snakemake@params[["conditions"]],
     plot_out = snakemake@output[["plot"]],
     stats_out = snakemake@output[["stats"]])

