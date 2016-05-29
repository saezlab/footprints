library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')

OUTFILE = commandArgs(TRUE)[1]

# load sanger data
scores = gdsc$mutated_genes(intogen=TRUE, drop=TRUE) + 0

save(scores, file=OUTFILE)
