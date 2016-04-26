library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
surv = import('./util')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/pathways_mapped/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "cont_speed_matrix.RData"

# load scores, only select primary tumors & map to patient IDs
scores = surv$load(file=INFILE)

pan_cov = surv$pancan(scores)
tissue = surv$tissue(scores)

save(pan_cov, tissue, file=OUTFILE)
