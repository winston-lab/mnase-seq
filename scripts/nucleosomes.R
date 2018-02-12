library(tidyverse)
library(forcats)
library(viridis)
library(ggthemes)

theme_default = theme_light() +
    theme(text = element_text(size=12, color="black", face="bold"),
          axis.text = element_text(size=12, color="black"),
          plot.title = element_text(size=12),
          plot.subtitle = element_text(size=10, color="black", face="plain"),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size=12, color="black", face="bold"),
          axis.title.y = element_text(angle=0, hjust=1, vjust=0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=12))

main = function(condition_label, control_label, in_table, shift_hist_out, occupancy_volcano_out,
                fuzziness_volcano_out, occupancy_v_fuzziness_out, occupancy_violin_out, occupancy_freqpoly_out,
                occupancy_ecdf_out, fuzziness_violin_out, fuzziness_freqpoly_out, fuzziness_ecdf_out){
    df = read_tsv(in_table, col_types = 'ciiiiiiiddddddddddddddd')
    
    #histogram of distance between nucleosome dyads (condition-control)
    shift_hist = ggplot(data = df, aes(x=abs(treat_smt_loca-control_smt_loca))) +
        geom_histogram(binwidth=10, boundary=0, fill="#114477", color="black") +
        xlab("nucleosome dyad shift (bp)") +
        ggtitle(paste("nucleosome dyad shifts,", condition_label, "vs.", control_label)) +
        theme_default +
        theme(axis.title.y = element_blank())
    ggsave(shift_hist_out, plot=shift_hist, width=16, height=8, units="cm")
    
    #occupancy volcano plot
    occupancy_volcano = ggplot(data = df, aes(x=point_log2FC, y=-log10(point_diff_FDR))) +
        stat_bin_hex(geom="point", aes(color=log10(..count..), fill=log10(..count..)),
                     shape=16, stroke=0, binwidth=c(0.02, 0.02), alpha=0.6) +
        scale_color_viridis(option="inferno", guide=FALSE) + 
        scale_fill_viridis(option="inferno", guide=FALSE) +
        xlab(bquote(bold(log[2] ~ frac(occupancy ~ .(condition_label),
                                       occupancy ~ .(control_label))))) +
        ylab(expression(bold(-log[10] ~ FDR))) +
        ggtitle(paste("nucleosome occupancy:", condition_label, "vs.", control_label)) +
        theme_default
    ggsave(occupancy_volcano_out, plot = occupancy_volcano, width=16, height=12, units="cm")
    
    #fuzziness volcano plot
    fuzziness_volcano = ggplot(data = df, aes(x=fuzziness_log2FC, y=-log10(fuzziness_diff_FDR))) +
        stat_bin_hex(geom="point", aes(color=log10(..count..), fill=log10(..count..)),
                     shape=16, stroke=0, binwidth=c(0.02, 0.02), alpha=0.6) +
        scale_color_viridis(option="inferno", guide=FALSE) + 
        scale_fill_viridis(option="inferno", guide=FALSE) +
        xlab(bquote(bold(log[2] ~ frac(fuzziness ~ .(condition_label),
                                       fuzziness ~ .(control_label))))) +
        ylab(expression(bold(-log[10] ~ FDR))) +
        ggtitle(paste("nucleosome fuzziness:", condition_label, "vs.", control_label),
                subtitle = expression(fuzziness %==% SD ~ of ~ distances ~ from ~ summit)) +
        theme_default
    ggsave(fuzziness_volcano_out, plot = fuzziness_volcano, width=16, height=12, units="cm")
    
    melted = df %>% select(control_smt_val, treat_smt_val, control_fuzziness_score, treat_fuzziness_score) %>% 
        gather(condition, occupancy, -c(control_fuzziness_score, treat_fuzziness_score)) %>% 
        mutate(condition = if_else(condition=="control_smt_val", control_label, condition_label)) %>% 
        mutate(condition=fct_inorder(condition, ordered=TRUE)) %>% 
        transmute(condition=condition, occupancy=occupancy,
                  fuzziness = if_else(condition==control_label, control_fuzziness_score, treat_fuzziness_score))
    
    #occupancy vs fuzziness
    occupancy_v_fuzziness = ggplot(data = melted, aes(x=occupancy, y=fuzziness)) +
        stat_bin_hex(geom="point", aes(color=log10(..count..)), binwidth=c(.1,.1), shape=16, stroke=0, alpha=0.6) +
        scale_x_log10() +
        scale_color_viridis(option="inferno", guide=FALSE) +
        scale_fill_viridis(option="inferno", guide=FALSE) +
        facet_wrap(~condition, dir="v") +
        ggtitle("nucleosome fuzziness vs. occupancy",
                subtitle = expression(fuzziness %==% SD ~ of ~ distances ~ from ~ summit)) +
        theme_default
    ggsave(occupancy_v_fuzziness_out, plot=occupancy_v_fuzziness, width=16, height=16, units="cm")
    
    #occupancy violin
    occupancy_violin = ggplot(data = melted, aes(x=condition, y=occupancy)) +
        geom_violin(bw=.01, draw_quantiles = .5, fill="#4477AA") +
        scale_y_log10() +
        ggtitle("nucleosome occupancy") +
        theme_default +
        theme(axis.title.x = element_blank())
    ggsave(occupancy_violin_out, plot=occupancy_violin, width=12, height=10, units="cm")
    
    #occupancy freqpoly
    occupancy_freqpoly = ggplot(data = melted, aes(x=occupancy, color=condition)) +
        geom_freqpoly(binwidth=.01, size=0.5) +
        scale_x_log10() +
        scale_color_ptol() +
        ggtitle("nucleosome occupancy") +
        theme_default +
        theme(axis.title.y = element_blank())
    ggsave(occupancy_freqpoly_out, plot=occupancy_freqpoly, width=16, height=8, units="cm")
    
    #occupancy eCDF
    occupancy_ecdf = ggplot(data = melted, aes(x=occupancy, color=condition)) +
        stat_ecdf(geom="step", size=0.5) +
        scale_x_log10() +
        scale_color_ptol() +
        ylab("cumulative probability") +
        ggtitle("eCDF of nucleosome occupancy") +
        theme_default +
        theme(axis.title.y = element_text(angle=90, hjust=0.5, size=10, face="plain"),
              legend.position = c(0.02,0.9),
              legend.justification = c(0,1),
              legend.background = element_blank(),
              legend.key = element_blank())
    ggsave(occupancy_ecdf_out, occupancy_ecdf, width=12, height=8, units="cm")
    
    #fuzziness violin
    fuzziness_violin = ggplot(data = melted, aes(x=condition, y=fuzziness)) +
        geom_violin(bw=.1, draw_quantiles = .5, fill="#4477AA") +
        ggtitle("nucleosome fuzziness",
                subtitle = expression(fuzziness %==% SD ~ of ~ distances ~ from ~ summit)) +
        theme_default +
        theme(axis.title.x = element_blank())
    ggsave(fuzziness_violin_out, plot=fuzziness_violin, width=12, height=10, units="cm")
    
    #fuzziness freqpoly
    fuzziness_freqpoly = ggplot(data = melted, aes(x=fuzziness, color=condition)) +
        geom_freqpoly(binwidth=.1, size=0.5) +
        scale_color_ptol() +
        ggtitle("nucleosome fuzziness",
                subtitle = expression(fuzziness %==% SD ~ of ~ distances ~ from ~ summit)) +
        theme_default +
        theme(axis.title.y = element_blank())
    ggsave(fuzziness_freqpoly_out, plot=fuzziness_freqpoly, width=16, height=8, units="cm")
    
    #fuzziness eCDF
    fuzziness_ecdf = ggplot(data = melted, aes(x=fuzziness, color=condition)) +
        stat_ecdf(geom="step", size=0.5) +
        scale_color_ptol() +
        ylab("cumulative probability") +
        ggtitle("eCDF of nucleosome fuzziness",
                subtitle = expression(fuzziness %==% SD ~ of ~ distances ~ from ~ summit)) +
        theme_default +
        theme(axis.title.y = element_text(angle=90, hjust=0.5, size=10, face="plain"),
              legend.position = c(0.02,0.9),
              legend.justification = c(0,1),
              legend.background = element_blank(),
              legend.key = element_blank())
    ggsave(fuzziness_ecdf_out, fuzziness_ecdf, width=12, height=8, units="cm")
}

main(condition_label = snakemake@wildcards[["condition"]],
     control_label = snakemake@wildcards[["control"]],
     in_table = snakemake@input[["results"]],
     shift_hist_out = snakemake@output[["shift_hist"]],
     occupancy_volcano_out = snakemake@output[["occupancy_volcano"]],
     fuzziness_volcano_out = snakemake@output[["fuzziness_volcano"]],
     occupancy_v_fuzziness_out = snakemake@output[["occupancy_v_fuzziness"]],
     occupancy_violin_out = snakemake@output[["occupancy_violin"]],
     occupancy_freqpoly_out = snakemake@output[["occupancy_freqpoly"]],
     occupancy_ecdf_out = snakemake@output[["occupancy_ecdf"]],
     fuzziness_violin_out = snakemake@output[["fuzziness_violin"]],
     fuzziness_freqpoly_out = snakemake@output[["fuzziness_freqpoly"]],
     fuzziness_ecdf_out = snakemake@output[["fuzziness_ecdf"]])
