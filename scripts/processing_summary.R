library(tidyverse)
library(forcats)
library(viridis)

survival_plot = function(df, scalefactor, ylabel){
    plot = ggplot(data = df, aes(x=step, y=count/scalefactor, group=sample)) +
        # geom_hline(aes(yintercept=count/scalefactor), color="grey50", size=0.2) +
        geom_step(direction="vh", position=position_nudge(x=0.5),
                  color="#114477", size=0.8) +
        scale_x_continuous(expand=c(0,0), breaks=2:4, name=NULL,
                           labels=c("raw reads", "reads cleaned",
                                    "aligned and paired")) +
        scale_y_continuous(sec.axis=dup_axis(), name=ylabel) +
        facet_grid(sample~., switch="y") +
        theme_light() +
        theme(strip.placement="outside",
              strip.background = element_blank(),
              text = element_text(size=12, color="black", face="bold"),
              strip.text.y = element_text(size=12, angle=-180, color="black", hjust=1),
              axis.text.x = element_text(size=12, color="black", face="bold", angle=30, hjust=0.95),
              axis.text.y = element_text(size=10, color="black", face="plain"),
              axis.title.y.right = element_blank(),
              panel.grid.major.x = element_line(color="grey40"),
              # panel.grid.major.y = element_blank(),
              # panel.grid.minor.y = element_blank(),
              plot.subtitle = element_text(size=12, face="plain"))
    return(plot)
}

main = function(in_table, surv_abs_out, surv_rel_out, loss_out){
    df = read_tsv(in_table) %>%
        mutate(sample=fct_inorder(sample, ordered=TRUE))
    
    nsamples = nrow(df)
    
    loss = df %>% gather(step, count, -sample, factor_key=TRUE) %>% 
        group_by(sample) %>% 
        mutate(og_count = lag(count)) %>% 
        filter(step != "raw") %>% 
        mutate(loss = (og_count-count)/og_count)
    
    #some hacking to get a survival-curve like thing
    #TODO: make the color fill the AUC?
    survival = df %>% mutate(dummy=raw) %>% 
        select(sample, dummy, 2:4) %>% 
        gather(step, count, -sample, factor_key=TRUE) %>% 
        mutate_at(vars(step), as.numeric)
    
    surv_abs = survival_plot(survival, scalefactor = 1e6, ylabel = "library size (M reads)") +
        ggtitle("read processing summary",
                subtitle = "absolute library size")
    
    surv_rel = survival_plot(survival %>% group_by(sample) %>% mutate(count=count/max(count)),
                             scalefactor = .01, ylabel = "% of raw reads") +
        ggtitle("read processing summary", subtitle = "relative to library size")
        
    ggsave(surv_abs_out, plot=surv_abs, width=14, height=2+2.5*nsamples, units="cm")
    ggsave(surv_rel_out, plot=surv_rel, width=14, height=2+2.5*nsamples, units="cm")
    
    loss_plot = ggplot(data = loss, aes(x=step, y=0, fill=loss)) +
        geom_raster() +
        geom_text(aes(label=round(loss, 2)), size=4) +
        scale_fill_viridis(name="% loss", guide=guide_colorbar(barheight = 10, barwidth=1)) +
        scale_color_viridis(guide=FALSE) +
        scale_x_discrete(labels = c("reads cleaned", "aligned",
                                    "aligned and paired"),
                         expand=c(0,0), name=NULL) +
        scale_y_continuous(breaks=0, expand=c(0,0), name=NULL) +
        facet_grid(sample~., switch="y") +
        ggtitle("read processing percent loss") +
        theme_light() +
        theme(strip.placement="outside",
              strip.background = element_blank(),
              text = element_text(size=12, color="black", face="bold"),
              strip.text.y = element_text(size=12, angle=-180, color="black", hjust=1),
              axis.text.x = element_text(size=12, color="black", face="bold", angle=30, hjust=0.95),
              axis.text.y = element_blank(),
              axis.title.y.right = element_blank(),
              plot.subtitle = element_text(size=12, face="plain"),
              panel.border = element_blank())
    
    ggsave(loss_out, plot=loss_plot, width=14, height=2+1.5*nsamples, units="cm")
}

main(in_table = snakemake@input[[1]],
     surv_abs_out = snakemake@output[["surv_abs_out"]],
     surv_rel_out = snakemake@output[["surv_rel_out"]],
     loss_out = snakemake@output[["loss_out"]])
