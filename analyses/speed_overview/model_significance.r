library(dplyr)
library(cowplot)
io = import('io')

top100 = function() {
    speed = io$load(module_file('../../model/model_matrix.RData'))$assocs %>%
        group_by(pathway) %>%
        top_n(100, -p.value)
}

plot = function(top) {
    top %>%
        group_by(pathway) %>%
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

if (is.null(module_name())) {
    pdf("model_significance.pdf")
    print(plot(top100()))
    dev.off()
}
