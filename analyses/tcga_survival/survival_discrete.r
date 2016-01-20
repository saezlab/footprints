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

# load and select primary tumors only, one sample per patient
scores = surv$load(file=INFILE)
clinical = surv$clinical

# make sure we have pathway scores and clinical on the same subset
ar$intersect(scores, clinical$barcode, along=1)
nnas = colSums(!is.na(scores))
if (any(nnas < 10)) {
    warning("Dropping with less than 10 values: ",
            paste(colnames(scores)[nnas >= 10], collapse=", "))
    scores = scores[,nnas >= 10]
}

# discretize in quartiles ("down", "unknown"*2, "up")
scores = scores %>%
    ar$map(along=1, subsets=clinical$study, surv$discretize_quartiles)

# calculate survival association using cox proportional hazards model
assocs.pan = surv$pancan(scores, clinical)
assocs.tissue = surv$tissue(scores, clinical)

# save the volcano plots in pdf
pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

assocs.pan %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1) %>%
    print()

fits = assocs.pan %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(5)
if (nrow(fits) >= 1)
    apply(fits, 1, surv$row2survFit)

assocs.tissue %>%
    plt$color$p_effect("adj.p", dir=-1, thresh=0.1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) %>%
    print()

fits = assocs.tissue %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(15)
if (nrow(fits >= 1))
    apply(fits, 1, surv$row2survFit)
