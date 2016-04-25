# get pathway scores and mutations, and correlate them with each other
library(dplyr)
b = import('base')
io = import('io')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "assocs_driver_mapped/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "plot_driver_mapped/speed_matrix.pdf"

assocs = io$load(INFILE)
subsets = c("pan_cov", "pan", setdiff(sort(unique(assocs$subset)), c("pan_cov", "pan")))

pdf(OUTFILE, paper="a4r", width=26, height=20)

for (subs in subsets) {
    if (grepl("pan", subs))
        pt_size = 0.5
    else
        pt_size = 5

    assocs %>%
        filter(subset == subs) %>%
        plt$color$p_effect(pvalue="adj.p", thresh=0.1) %>%
        plt$volcano(base.size=pt_size, p=0.1) +
            ggtitle(subs)
}

dev.off()
