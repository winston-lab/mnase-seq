library(tidyverse)
library(GGally)
library(viridis)

main = function(intable, binsize, pcount, samplelist, outpath){
    df = intable %>% read_tsv() %>%
        gather(key=sample, value=signal, -name) %>%
        filter(sample %in% samplelist) %>%
        mutate_at(vars(sample), ~(fct_inorder(., ordered=TRUE))) %>%
        spread(sample, signal) %>%
        select(-name)

    df = df[which(rowSums(df)>0),]

    cor_matrix = df %>% na_if(0) %>% log10() %>%
        cor(method="pearson", use="pairwise.complete.obs")

    maxsignal = max(df) + pcount
    mincor = min(cor_matrix) * 0.98
    plots = list()

    #for each row
    for (i in 1:ncol(df)){
        #for each column
        for (j in 1:ncol(df)){
            idx = ncol(df)*(i-1)+j
            if (i < j){
                #upper right (correlation)
                cor_value = cor_matrix[i,j]
                plot = ggplot(data = tibble(x=c(0,1), y=c(0,1), corr=cor_value)) +
                        geom_rect(aes(fill=corr), xmin=0, ymin=0, xmax=1, ymax=1) +
                        annotate("text", x=0.5, y=0.5, label=sprintf("%.2f",round(cor_value,2)), size=10*abs(cor_value)) +
                        scale_x_continuous(breaks=NULL) +
                        scale_y_continuous(breaks=NULL) +
                        scale_fill_distiller(palette="Blues", limits = c(mincor,1), direction=1)
                plots[[idx]] = plot
            } else if (i == j){
                #top left to bot right diag (density)
                subdf = df %>% select(i) %>% gather(sample, value)
                plot = ggplot(data = subdf, aes(x=(value+pcount))) +
                        geom_density(aes(y=..scaled..), fill="#114477", size=0.8) +
                        scale_y_continuous(breaks=c(0,.5,1)) +
                        scale_x_log10(limit = c(pcount, maxsignal)) +
                        annotate("text", x=.90*maxsignal, y=0.5, hjust=1,
                                 label=unique(subdf$sample), size=4, fontface="bold")
                plots[[idx]] = plot
            } else {
                #bottom left (scatter)
                #filtering is an optional hack to avoid the (0,0) bin taking up
                #all of the colorspace
                subdf = df %>% select(i,j) %>% gather(xsample, xvalue, -1) %>%
                            gather(ysample, yvalue, -c(2:3)) #%>%
                            # filter(!(xvalue < 6*pcount & yvalue < 6*pcount))
                plot = ggplot(data = subdf, aes(x=xvalue+pcount, y=yvalue+pcount)) +
                            geom_abline(intercept = 0, slope=1, color="grey80", size=.5) +
                            stat_bin_hex(geom="point", aes(color=log10(..count..)), binwidth=c(.04,.04), size=.5, shape=16, stroke=0) +
                            scale_fill_viridis(option="inferno") +
                            scale_color_viridis(option="inferno") +
                            scale_x_log10(limit = c(pcount, maxsignal)) +
                            scale_y_log10(limit = c(pcount, maxsignal))
                plots[[idx]] = plot
            }
        }
    }

    mat = ggmatrix(plots, nrow=ncol(df), ncol=ncol(df),
                   title = paste0("MNase-seq dyad signal, ", binsize, "bp bins"),
                   xAxisLabels = names(df), yAxisLabels = names(df), switch="both") +
                    theme_light() +
                    theme(plot.title = element_text(size=12, color="black", face="bold"),
                          axis.text = element_text(size=9),
                          strip.background = element_blank(),
                          strip.text = element_text(size=12, color="black", face="bold"),
                          strip.text.y = element_text(angle=180, hjust=1),
                          strip.placement="outside",
                          strip.switch.pad.grid = unit(0, "points"),
                          strip.switch.pad.wrap = unit(0, "points"))
    w = 3+ncol(df)*4
    h = 9/16*w+0.5
    ggsave(outpath, mat, width=w, height=h, units="cm")
    print(warnings())
}

main(intable = snakemake@input[[1]],
     binsize = snakemake@wildcards[["windowsize"]],
     pcount = snakemake@params[["pcount"]],
     samplelist = snakemake@params[["samplelist"]],
     outpath = snakemake@output[[1]])

