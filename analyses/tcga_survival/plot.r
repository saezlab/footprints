library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
surv = import('./util')

INFILE = commandArgs(TRUE)[1] %or% "assocs_cont_mapped/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "plots_cont_mapped/speed_matrix.pdf"

# load scores, only select primary tumors & map to patient IDs
assocs = io$load(file=INFILE)

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

assocs$pan_cov %>%
    plt$color$p_effect("adj.p", dir=-1, thresh=0.1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1, p=0.1) %>%
    print()

assocs$tissue %>%
    plt$color$p_effect("adj.p", dir=-1, thresh=0.1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1, label_top=50) %>% #+ ggtitle(sum(clinical$adj.p < 0.1)) %>%
    print()
