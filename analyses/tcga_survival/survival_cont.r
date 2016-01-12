library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
util = import('./util')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_norm.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "cont_speed_norm.pdf"

scores = io$load(INFILE)

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

util$pancan(scores) %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1) %>%
    print()

util$tissue(scores) %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) + ggtitle(sum(clinical$adj.p < 0.1)) %>%
    print()
