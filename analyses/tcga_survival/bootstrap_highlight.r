library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
util = import('./util')

bootstrap = function(tissue, pathway, scores, clinical, n=100) {
    one_sample = function(seed=12507) {
        s = scores %>%
            tcga$filter(tissue=tissue)
        s = scores[sample(rownames(s), replace=TRUE),pathway,drop=FALSE]
        rownames(s) = substr(rownames(s), 1, 12)
        cc = clinical[match(rownames(s), clinical$barcode),]
        rownames(s) = make.names(rownames(s), unique=TRUE)
        cc$barcode = make.names(cc$barcode, unique=TRUE)
        util$tissue(s, cc)
    }

    lapply(12507 + 1:100, one_sample) %>%
        bind_rows()
}

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/pathways_mapped/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "bootstrap_highlight.RData"

# load and select primary tumors only, one sample per patient
scores = io$load(file=INFILE)
clinical = util$clinical

# make sure we have pathway scores and clinical on the same subset
tcga$intersect(scores, clinical$barcode, along=1)
nnas = colSums(!is.na(scores))
if (any(nnas < 10)) {
    warning("Dropping with less than 10 values: ",
            paste(colnames(scores)[nnas >= 10], collapse=", "))
    scores = scores[,nnas >= 10]
}

# discretize in quartiles ("down", "unknown"*2, "up")
scores = scores %>%
    ar$map(along=1, subsets=clinical$study, util$discretize_quartiles)

bootstrap = bind_rows(list(
    bootstrap("KIRC", "TNFa", scores, clinical, 100),
    bootstrap("LGG", "JAK-STAT", scores, clinical, 100),
    bootstrap("ACC", "p53", scores, clinical, 100)
))

save(bootstrap, file=OUTFILE)
