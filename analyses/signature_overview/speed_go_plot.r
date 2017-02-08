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
plot_hypergeom = function(mat) {
}

if (is.null(module_name())) {
    INFILE = commandArgs(TRUE)[1] %or% "speed_go_hypergeom.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "speed_go_plot.pdf"

    mat = io$load(INFILE)

    pdf(file=OUTFILE)
    on.exit(dev.off)
    plot_hypergeom(mat)
}
