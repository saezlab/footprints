library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
util = import('./util')

OUTFILE = commandArgs(TRUE)[2] %or% "pancan_heatmap.pdf"

fnames = c("speed_matrix", "gsea_reactome", "gsea_go", "spia", "pathifier")
assocs = lapply(fnames, function(fn) {
    fname = paste0("../../scores/tcga/", fn, ".RData")
    scores = io$load(fname)
    util$pancan(fname, util$clinical) %>%
}) %>%
    df$add_col_names("method") %>%
    bind_rows() %>%
    mutate(estimate = ifelse(adj.p < 0.1, estimate, NA)) %>%
    mutate(label = ifelse(adj.p < 0.05, "*", "")) %>%
    mutate(label = ifelse(adj.p < 0.01, "***", label))

pdf(OUTFILE, paper="a4r", width=10, height=5)
on.exit(dev.off)

plt$matrix(assocs, estimate ~ method + scores) +
    xlab("Pathway") +
    ylab("Method") +
    ggtitle("Pan-cancer survival associations")
