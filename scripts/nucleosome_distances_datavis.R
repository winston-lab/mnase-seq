library(tidyverse)
library(ggthemes)
library(ggridges)

theme_default= theme_light() +
    theme(axis.text=element_text(color="black"),
          axis.title.x.top=element_blank(),
          panel.grid.minor.x=element_blank(),
          strip.background=element_blank(),
          strip.text=element_text(color="black",
                                  hjust=0))

import = function(tsv_path,
                  condition_id){
    read_tsv(tsv_path) %>%
        mutate(condition=condition_id) %>%
        return()
    }

main = function(
    control_tsv_path = "control.tsv",
    condition_tsv_path = "condition.tsv",
    control_id = "WT-37C",
    condition_id = "spt6-YW-37C",
    ridgelines_out="ridgelines.pdf",
    jitter_out="jitter.pdf"){

    df = bind_rows(import(control_tsv_path, control_id),
                   import(condition_tsv_path, condition_id)) %>%
        mutate(nuc_index_jitter = jitter(nuc_index,amount=0.4),
               condition=ordered(condition, levels=c(control_id, condition_id)))

    summary_df = df %>%
        group_by(condition, nuc_index) %>%
        summarise(distance_mean=mean(distance_from_previous),
                  distance_sd=sd(distance_from_previous))

    ridgelines = ggplot(data=filter(df, nuc_index<10),
          aes(x=distance_from_previous,
              y=nuc_index,
              group=interaction(nuc_index, condition),
              fill=condition,
              color=condition)) +
       geom_density_ridges2(alpha=0.2,
                            size=0.3,
                            rel_min_height=0.01,
                           scale=0.9,
                           quantile_lines=TRUE,
                           quantiles=2,
                           panel_scaling = FALSE,
                           bandwidth=3) +
       scale_y_discrete(limits=seq(-0.5,9.5, 0.5),
                          breaks=seq(-0.5,9.5),
                          labels=c("TSS", paste("+", seq(1, 10))),
                          name="nucleosomes") +
       scale_x_continuous(limits=c(NA, 332),
                          name="distance (bp)",
                          breaks=scales::pretty_breaks(6)) +
       scale_fill_ptol(guide=guide_legend(label.position="top",
                                          keyheight=unit(10, "pt"))) +
       scale_color_ptol() +
       theme_default +
       theme(panel.grid.major.y=element_blank(),
             panel.grid.major.x=element_line(size=0.1),
             legend.position="top",
             legend.title=element_blank(),
             legend.spacing.y=unit(2, "pt"),
             legend.margin=margin(b=-8, unit="pt"))

    ggsave(ridgelines_out,
          plot=ridgelines,
          width=16,
          height=9,
          units="cm")

    jitter = ggplot(data=df,
          aes(x=nuc_index_jitter,
              y=distance_from_previous)) +
       stat_binhex(geom="point",
                   aes(color=..count..),
                   binwidth=c(0.08,2.49),
                   alpha=0.65,
                   shape=16,
                   size=0.6) +
       geom_text(data=summary_df,
                 aes(x=nuc_index,
                     label=paste("textstyle(atop(phantom(.) * ", round(distance_mean), ", phantom(.) %+-%", round(distance_sd), "))"),
                     y=-20),
                 parse=TRUE,
                 size=10/72*25.4,
                 vjust=0.5) +
       facet_grid(.~condition) +
       scale_x_continuous(limits=c(NA, 9.5),
                          breaks=seq(-0.5,9.5),
                          labels=c("TSS", paste("+", seq(1, 10))),
                          name="nucleosomes") +
       scale_y_continuous(limits=c(NA, 332),
                          name="distance (bp)") +
       scale_color_viridis_c() +
       theme_light() +
       theme(legend.position="none",
             axis.text=element_text(color="black"),
             axis.title.x.top=element_blank(),
             panel.grid.minor.x=element_blank(),
             strip.background=element_blank(),
             strip.text=element_text(color="black",
                                     hjust=0))

    ggsave(jitter_out,
          plot=jitter,
          width=16,
          height=9,
          units="cm")
}

main(control_tsv_path = snakemake@input[["control"]],
     condition_tsv_path = snakemake@input[["condition"]],
     control_id = snakemake@wildcards[["control"]],
     condition_id = snakemake@wildcards[["condition"]],
     ridgelines_out= snakemake@output[["ridgelines"]],
     jitter_out= snakemake@output[["jitter"]])

# violins = ggplot(data=df,
#        aes(x=nuc_index,
#            y=distance_from_previous,
#            group=interaction(nuc_index, condition),
#            color=condition,
#            fill=condition)) +
#     geom_violin(position=position_dodge(width=0.8),
#                 bw=10,
#                 size=0.2,
#                 alpha=0.9,
#                 width=1) +
#     geom_boxplot(outlier.shape=NA,
#                  notch=TRUE,
#                  width=0.15,
#                  fill="black",
#                  color="black",
#                  size=0.1,
#                  alpha=0.2,
#                  position=position_dodge(width=0.8)) +
#     scale_x_continuous(limits=c(NA, 9.5),
#                        breaks=seq(-0.5,9.5),
#                        labels=c("TSS", paste("+", seq(1, 10))),
#                        name="nucleosomes",
#                        expand=c(0,0.2)) +
#     scale_y_continuous(limits=c(NA, 332),
#                        name="distance (bp)",
#                        breaks=scales::pretty_breaks(6)) +
#     scale_color_ptol(guide=guide_legend(label.position="top",
#                                         keyheight=unit(8, "pt"))) +
#     scale_fill_ptol() +
#     theme_default +
#     theme(legend.position="top",
#           legend.title=element_blank(),
#           legend.spacing.y=unit(2,"pt"),
#           legend.margin=margin(b=-10, unit="pt"),
#           panel.grid.major.x=element_line(linetype="dashed"),
#           panel.grid.minor.y=element_blank(),
#           panel.grid.major.y=element_blank())
#
# ggsave("violins.pdf",
#        plot=violins,
#        width=16,
#        height=9,
#        units="cm")

# ggplot(data=df,
#        aes(x=nuc_index,
#            y=distance_from_previous,
#            group=transcript_name)) +
#     geom_line(alpha=0.01) +
#     scale_x_continuous(limits=c(NA, 11),
#                        breaks=seq(-0.5,10.5),
#                        labels=c("TSS", paste("+", seq(1, 11)))) +
#     scale_y_continuous(limits=c(NA, 500))

# dd = df %>%
#     mutate(condition=if_else(condition==condition_id, "condition", "control")) %>%
#     pivot_wider(id_cols=c(chrom, transcript_start, transcript_end, transcript_name, strand, nuc_index),
#                 names_from=condition,
#                 values_from=distance_from_TSS)
#
# ggplot(data=filter(dd, nuc_index<10),
#        aes(x=control,
#            y=condition,
#            group=transcript_name)) +
#     geom_abline(slope=1,
#                 intercept=0,
#                 color="black") +
#     geom_point(aes(color=factor(nuc_index)),
#                size=0.5,
#                shape=1,
#                alpha=0.1) +
#     # geom_line(alpha=0.01) +
#     coord_fixed()
