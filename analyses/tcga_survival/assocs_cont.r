library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
surv = import('./util')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "cont_speed_matrix.RData"

# load scores, only select primary tumors & map to patient IDs
scores = surv$load(file=INFILE)

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

pancan = surv$pancan(scores)
tissue = surv$tissue(scores)

save(pancan, tissue, file=OUTFILE)
