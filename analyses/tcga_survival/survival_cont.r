library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
surv = import('./util')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "cont_speed_matrix.pdf"

# load scores, only select primary tumors & map to patient IDs
scores = surv$load(file=INFILE)

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

surv$pancan(scores) %>%
    plt$color$p_effect("adj.p", dir=-1, thresh=0.1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1, p=0.1) %>%
    print()

surv$tissue(scores) %>%
    plt$color$p_effect("adj.p", dir=-1, thresh=0.1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) %>% #+ ggtitle(sum(clinical$adj.p < 0.1)) %>%
    print()
