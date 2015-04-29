#!/usr/bin/env Rscript
b = import('base')
io = import('io')
ma = import('process/microarray')

ACCESSION = commandArgs(TRUE)[1] %or% 'E-GEOD-8346'
OUTFILE = commandArgs(TRUE)[2] %or% paste0(ACCESSION, ".RData")

# read raw data, normalize, qc
expr = ArrayExpress::ArrayExpress(ACCESSION) %>%
    ma$qc() %>%
    ma$normalize()

# save result
save(expr, file=OUTFILE)
