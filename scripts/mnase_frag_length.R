library(tidyverse)
library(forcats)

main = function(in_table, out_path){
    df = read_tsv(in_table) %>%
        gather(key=sample, value=count, -fragsize) %>%
        mutate(sample = fct_inorder(sample, ordered=TRUE)) %>%
        group_by(sample) %>%
        mutate(density = count/sum(count))

    plot = ggplot(data = df, aes(x=fragsize, y=density)) +
        geom_area(fill="#114477", color="black") +
        facet_grid(sample~., switch="y") +
        scale_y_continuous(breaks = scales::pretty_breaks(n=2)) +
        xlab("fragment size (bp)") +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(color="black"),
              axis.text.x = element_text(size=12),
              axis.text.y = element_text(face="plain"),
              strip.background = element_blank(),
              strip.text = element_text(color="black", size=12),
              strip.placement = "outside",
              strip.text.y = element_text(angle=-180, hjust=1))

    ggsave(out_path, plot=plot,
           width=14, height=2+1.5*n_distinct(df[["sample"]]),
           units="cm", limitsize=FALSE)
}

main(in_table = snakemake@input[["table"]],
     out_path = snakemake@output[["plot"]])
