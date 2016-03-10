library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
surv = import('./util')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_matrix.RData"

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
pan = surv$pancan(scores, clinical)
tissue = surv$tissue(scores, clinical)

save(pan, tissue, file=OUTFILE)
