library(dplyr)
io = import('io')
plt = import('plot')

plot_pval = function() {
    io$load(module_file('../../model/model_matrix.RData'))$assocs %>%
        group_by(pathway) %>%
        top_n(100, -p.value) %>%
        mutate(i=match(1:length(p.value), order(p.value))) %>%
        ungroup() %>%
        ggplot(aes(x=i, y=-log10(adj.p))) +
        geom_area(color="blue", fill="blue", alpha=0.5) +
        facet_wrap(~pathway) +
        xlab("Top 100 genes") +
        ylab("- log(FDR)") +
        scale_y_continuous(trans="log10") +
        geom_hline(yintercept=-log10(0.05), lwd=0.3, linetype='dashed') +
        geom_hline(yintercept=-log10(1e-10), lwd=0.3, linetype='dotted')
}

plot_zscore = function() {
    io$load(module_file('../../model/model_matrix.RData'))$assocs %>%
        group_by(pathway) %>%
        top_n(100, -p.value) %>%
        mutate(i=match(1:length(zscore), order(-zscore))) %>%
        ungroup() %>%
        ggplot(aes(x=i, y=zscore)) +
        geom_area(color="red", fill="red", alpha=0.5) +
        facet_wrap(~pathway) +
        xlab("Top 100 genes") +
        ylab("z coefficient")
}

if (is.null(module_name())) {
    OUTFILE = commandArgs(TRUE)[1]

    pdf(OUTFILE)
    print(plot_pval())
    print(plot_zscore())
    dev.off()
}
