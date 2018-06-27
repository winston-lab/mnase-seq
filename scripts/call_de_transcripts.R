library(tidyverse)
library(forcats)
library(DESeq2)
library(viridis)
library(scales)
library(gridExtra)
library(ggrepel)

get_countdata = function(path, samples){
    df = read_tsv(path, col_names=TRUE) %>% select(c("name", samples)) %>%
            column_to_rownames(var="name") %>% as.data.frame()
    df = df[rowSums(df)>1,]
    return(df)
}

mean_sd_plot = function(df, ymax){
    ggplot(data = df, aes(x=rank, y=sd)) +
        geom_hex(aes(fill=log10(..count..), color=log10(..count..)), bins=40, size=0) +
        geom_smooth(color="#4292c6") +
        scale_fill_viridis(option="inferno", name=expression(log[10](count)), guide=FALSE) +
        scale_color_viridis(option="inferno", guide=FALSE) +
        scale_x_continuous(trans=reverselog_trans(10), name="rank(mean expression)") +
        scale_y_continuous(limits = c(NA, ymax)) +
        ylab("SD") +
        theme_light() +
        theme(text = element_text(size=8))
}

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

separate_name = function(df){
    df %>% separate(name, into=c('transcript_id','chrstrand','start','end'), sep="~", convert=TRUE) %>%
        separate(chrstrand, into=c('chrom', 'strand'), sep="-") %>%
        mutate_at(vars(strand), funs(if_else(.=="minus", "-", "+"))) %>%
        return()
}

call_de_bases = function(intable, norm, sitable, samples, groups, condition, control, alpha, lfc,
                         results_all, results_up, results_down, results_unch,
                         bed_all, bed_up, bed_down, bed_unch,
                         normcounts, rldcounts, qcplots){
    #import data
    countdata = get_countdata(intable, samples)
    coldata = data.frame(condition=factor(groups,
                                          levels = c(control, condition)),
                         row.names=names(countdata))
    #run DESeq2
    dds = DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ condition)

    if (norm=="spikenorm"){
        sicountdata = get_countdata(sitable, samples)
        sidds = DESeqDataSetFromMatrix(countData = sicountdata,
                                       colData = coldata,
                                       design = ~ condition)
        sidds = sidds %>% estimateSizeFactors()
        sizeFactors(dds) = sizeFactors(sidds)
    } else {
        dds = dds %>% estimateSizeFactors()
    }
    dds = dds %>% estimateDispersions() %>% nbinomWaldTest()

    #extract normalized counts and write to file
    ncounts = dds %>% counts(normalized=TRUE) %>% as.data.frame() %>%
        rownames_to_column(var='name') %>%
        as_tibble() %>%
        separate_name() %>%
        mutate_if(is.numeric, round, 3) %>%
        write_tsv(path=normcounts, col_names=TRUE)

    ncountsavg = ncounts %>%
        gather(sample, value, -c(transcript_id, chrom, strand, start, end)) %>%
        mutate(group = if_else(sample %in% samples[groups==condition], condition, control)) %>%
        group_by(transcript_id, chrom, strand, start, end, group) %>%
        summarise(mean = mean(value)) %>% spread(group, mean) %>%
        ungroup()

    #plot sd vs. mean for unshrunken (log2) counts
    ntd = dds %>% normTransform() %>% assay() %>% as.data.frame() %>%
        rownames_to_column(var="name") %>% as_tibble() %>%
        gather(sample, signal, -name) %>% group_by(name) %>%
        summarise(mean=mean(signal), sd=sd(signal)) %>%
        mutate(rank = min_rank(desc(mean)))
    maxsd = max(ntd$sd)*1.01
    ntdplot = mean_sd_plot(ntd, maxsd) +
        ggtitle(expression(paste("raw ", log[2]("counts"))))

    #extract rlog transformed counts and write to file
    rlogcounts = dds %>% rlog(blind=FALSE) %>% assay() %>% as.data.frame() %>%
        rownames_to_column(var="name") %>% as_tibble() %>%
        separate_name() %>%
        mutate_if(is.numeric, round, 3) %>%
        write_tsv(path=rldcounts, col_names=TRUE)

    #plot sd vs. mean for rlog transformed counts
    rld = rlogcounts %>%
        gather(sample, signal, -c(transcript_id, chrom, strand, start, end)) %>%
        group_by(transcript_id, chrom, strand, start, end) %>%
        summarise(mean=mean(signal), sd=sd(signal)) %>%
        ungroup() %>%
        mutate(rank = min_rank(desc(mean)))
    rldplot = mean_sd_plot(rld, maxsd) +
                ggtitle(expression(paste("regularized ", log[2]("counts"))))

    #extract DESeq2 results and write to file
    resdf = results(dds, alpha=alpha, lfcThreshold=lfc, altHypothesis="greaterAbs") %>%
        as_tibble() %>%
        rownames_to_column(var='name') %>%
        separate_name() %>%
        arrange(padj) %>%
        inner_join(ncountsavg, by=c('transcript_id', 'chrom', 'strand', 'start', 'end')) %>%
        mutate_at(c('pvalue','padj'), funs(-log10(.))) %>%
        mutate_if(is.numeric, round, 3) %>% dplyr::rename(logpval=pvalue, logpadj=padj, meanExpr=baseMean) %>%
        mutate_at(vars(start, end), funs(as.integer(.))) %>%
        write_tsv(path=results_all, col_names=TRUE)

    #plot library size vs sizefactor
    sfdf = dds %>% sizeFactors() %>% as_tibble() %>%
        rownames_to_column(var="sample") %>% dplyr::rename(sizefactor=value) %>%
        inner_join(colSums(countdata) %>% as_tibble() %>%
                       rownames_to_column(var="sample") %>% dplyr::rename(libsize=value),
                   by="sample")
    sfplot = ggplot(data = sfdf, aes(x=libsize/1e6, y=sizefactor)) +
        geom_smooth(method="lm", se=FALSE, color="red", size=0.5) +
        geom_point(size=0.5) +
        geom_text_repel(aes(label=sample), size=2) +
        xlab("library size (M reads)") +
        ylab("size factor (median of ratios)") +
        theme_light() +
        theme(text = element_text(size=8))

    #MA plot for differential expression
    resdf.sig = resdf %>% filter(logpadj> -log10(alpha))
    resdf.nonsig = resdf %>% filter(logpadj<= -log10(alpha))
    maplot = ggplot() +
        geom_hline(yintercept = 0, color="black", linetype="dashed") +
        geom_point(data = resdf.nonsig, aes(x=meanExpr, y=log2FoldChange),
                   color="black", alpha=0.3, stroke=0, size=0.7) +
        geom_point(data = resdf.sig, aes(x=meanExpr, y=log2FoldChange),
                   color="red", alpha=0.3, stroke=0, size=0.7) +
        scale_x_log10(name="mean of normalized counts") +
        ylab(substitute(log[2]~frac(cond,cont), list(cond=condition, cont=control))) +
        theme_light() +
        theme(text = element_text(size=8))

    volcano = ggplot() +
        geom_point(data = resdf.nonsig, aes(x=log2FoldChange, y = logpadj),
                   alpha=0.3, stroke=0, size=0.7) +
        geom_point(data = resdf.sig, aes(x=log2FoldChange, y = logpadj),
                   alpha=0.3, stroke=0, size=0.7) +
        geom_hline(yintercept = -log10(alpha), color="red", linetype="dashed") +
        xlab(substitute(log[2]~frac(cond,cont), list(cond=condition, cont=control))) +
        ylab(expression(-log[10]("p value"))) +
        theme_light() +
        theme(text = element_text(size=8))

    out = grid.arrange(sfplot, ggplot()+theme_void(),
                       ntdplot, rldplot,
                       maplot, volcano, ncol=2,
                       heights = unit(c(4, 6, 6), rep("cm",3)))
    ggsave(qcplots, out, height=18, width = 16, units="cm")


    #write out results and bed files for transcripts going up and down
    resdf %>% unite(col="score", log2FoldChange, logpadj, sep=":") %>%
        select(chrom, start, end, transcript_id, score, strand) %>%
        write_tsv(bed_all, col_names=FALSE)

    resdf %>% filter(logpadj > -log10(alpha) & log2FoldChange > 0) %>%
        write_tsv(results_up) %>%
        unite(col="score", log2FoldChange, logpadj, sep=":") %>%
        select(chrom, start, end, transcript_id, score, strand) %>%
        write_tsv(bed_up, col_names=FALSE)

    resdf %>% filter(logpadj > -log10(alpha) & log2FoldChange < 0) %>%
        write_tsv(results_down) %>%
        unite(col="score", log2FoldChange, logpadj, sep=":") %>%
        select(chrom, start, end, transcript_id, score, strand) %>%
        write_tsv(bed_down, col_names=FALSE)

    resdf %>% filter(logpadj <= -log10(alpha)) %>%
        write_tsv(results_unch) %>%
        unite(col="score", log2FoldChange, logpadj, sep=":") %>%
        select(chrom, start, end, transcript_id, score, strand) %>%
        write_tsv(bed_unch, col_names=FALSE)
}

qc = call_de_bases(intable = snakemake@input[["expcounts"]],
                   norm = snakemake@wildcards[["norm"]],
                   sitable = snakemake@input[["sicounts"]],
                   samples = snakemake@params[["samples"]],
                   groups = snakemake@params[["groups"]],
                   condition = snakemake@wildcards[["condition"]],
                   control = snakemake@wildcards[["control"]],
                   alpha = snakemake@params[["alpha"]],
                   lfc = snakemake@params[["lfc"]],
                   results_all = snakemake@output[["results_all"]],
                   results_up = snakemake@output[["results_up"]],
                   results_down = snakemake@output[["results_down"]],
                   results_unch = snakemake@output[["results_unch"]],
                   bed_all = snakemake@output[["bed_all"]],
                   bed_up = snakemake@output[["bed_up"]],
                   bed_down = snakemake@output[["bed_down"]],
                   bed_unch = snakemake@output[["bed_unch"]],
                   normcounts = snakemake@output[["normcounts"]],
                   rldcounts = snakemake@output[["rldcounts"]],
                   qcplots = snakemake@output[["qcplots"]])
