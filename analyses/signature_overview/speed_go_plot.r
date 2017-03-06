library(tidyverse)
b = import('base')
io = import('io')

#' Plot enrichment for piano results
#'
#' @param index  A data.frame with columns: pathway, set, p[^A]
plot_piano = function(index) {
    index = io$load(INFILE) %>%
        tidyr::gather(type, p.value, -pathway, -set) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(p.value) %>%
        group_by(pathway, set, type) %>%
        top_n(5) %>%
        ungroup()
}

#' Plot enrichment for hypergeometric test
#'
#' @param mat  A matrix of p-values
plot_hypergeom = function() {
    mat = io$load(module_file("speed_go_hypergeom.RData", mustWork=TRUE))

    df = reshape2::melt(mat) %>%
        select(set=Var1, pathway=Var2, p.value=value) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               pathway = factor(pathway, levels=gtools::mixedsort(
                    as.character(unique(pathway))))) %>%
        arrange(p.value) %>%
        group_by(pathway) %>%
        top_n(10, -p.value) %>%
        ungroup() %>%
        mutate(set = reorder(set, -adj.p),
               adj.p = -log10(adj.p))

    ggplot(df, aes(x=set, y=adj.p)) +
        geom_bar(stat='identity', fill="lightblue") +
        geom_hline(yintercept=-log10(0.1), linetype="dotted") +
        coord_flip() +
        geom_text(aes(y=0, label=paste(" ", set)), hjust=0) +
        facet_wrap(~ pathway, scales="free", ncol=2) +
        theme_minimal() +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) +
        xlab("Gene Ontology category") +
        ylab("- log10 (FDR)")
}

if (is.null(module_name())) {
    INFILE = commandArgs(TRUE)[1] %or% "speed_go_hypergeom.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "speed_go_plot.pdf"

    pdf(file=OUTFILE)
    on.exit(dev.off)
    plot_hypergeom()
}
