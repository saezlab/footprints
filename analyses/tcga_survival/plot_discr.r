library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
surv = import('./util')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_matrix.pdf"

load(INFILE)

# save the volcano plots in pdf
pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

pan %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1) %>%
    print()

fits = pan %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(5)
if (nrow(fits) >= 1)
    apply(fits, 1, surv$row2survFit)

tissue %>%
    plt$color$p_effect("adj.p", dir=-1, thresh=0.1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) %>%
    print()

fits = tissue %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(15)
if (nrow(fits >= 1))
    apply(fits, 1, surv$row2survFit)
