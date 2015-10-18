library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_matrix.pdf"

assocfile2plot = function(fname) {
    assocs.tissue = io$load(fname)

    # volcano plot for tissue subsets
    assocs.tissue %>%
        mutate(label = paste(Yf, scores, sep=":")) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
        plt$volcano(p=0.2)
}

pdf(OUTFILE, paper="a4r", width=26, height=20)
print(assocfile2plot("speed_matrix.RData") + ggtitle("response genes"))
print(assocfile2plot("gsea_reactome.RData") + ggtitle("reactome expr"))
dev.off()
