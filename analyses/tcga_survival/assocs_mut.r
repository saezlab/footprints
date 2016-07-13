library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
util = import('./util')

INFILE = commandArgs(TRUE)[1] %or% "../tcga_pathway_per_mutation/mut_driver_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "assocs_discr_all/mutation.RData"

# load and select primary tumors only, one sample per patient
scores = io$load(file=INFILE)
clinical = util$clinical

# filter by min number of mutations
scores = scores[, colSums(scores) > 10]
rownames(scores) = substr(rownames(scores), 1, 12)
ar$intersect(scores, clinical$barcode, along=1)

# calculate survival association using cox proportional hazards model
pan_cov = util$pancan(scores, clinical)
tissue = util$tissue(scores, clinical)

save(pan_cov, tissue, file=OUTFILE)
