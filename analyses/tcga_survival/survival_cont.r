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

# load scores, only select primary tumors & map to patient IDs
scores = io$load(INFILE)
scores = scores[substr(rownames(scores), 14, 16) == "01A",]
rownames(scores) = substr(rownames(scores), 1, 12)

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

#    plt$color$p_effect("p.value", dir=-1) %>%
util$pancan(scores) %>%
    plt$color$p_effect("adj.p", dir=-1, thresh=0.1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1, p=0.1) %>%
    print()

#    plt$color$p_effect("p.value", dir=-1, thresh=0.1) %>%
util$tissue(scores) %>%
    plt$color$p_effect("adj.p", dir=-1, thresh=0.1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) %>% #+ ggtitle(sum(clinical$adj.p < 0.1)) %>%
    print()
