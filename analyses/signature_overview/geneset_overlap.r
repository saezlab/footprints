library(dplyr)
library(grid)
library(gridExtra)
b = import('base')
io = import('io')
ar = import('array')
plt = import('plot')
config = import('../../config/methods.yaml')

set_overlap_matrix = function() {
    # import all sets from methods listed in config

    # compute their overlap

    # return a ggplot upper diagonal matrix
}

if (is.null(module_name())) {
    pdf("geneset_overlap.pdf", width=19, height=15)
    on.exit(dev.off)

    set_overlap_matrix()
}
