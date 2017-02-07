library(dplyr)
library(magrittr)
plt = import('plot')
get_genesets = import('../../util/genesets')$get_genesets

go_enrichment = function() {
}

bar_plot = function() {
}

if (is.null(module_name())) {
    pdf("geneset_overlap.pdf")
    on.exit(dev.off)

    sets = get_genesets()
    geneset_overlap_matrix(sets)
}
