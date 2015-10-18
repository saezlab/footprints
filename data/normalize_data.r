#!/usr/bin/env Rscript
b = import('base')
io = import('io')
ma = import('process/microarray')
idmap = import('process/idmap')

ACCESSION = commandArgs(TRUE)[1] %or% 'E-GEOD-8346'
OUTFILE = commandArgs(TRUE)[2] %or% paste0(ACCESSION, ".RData")

# read raw data, normalize, qc
expr = ArrayExpress::ArrayExpress(ACCESSION, drop=FALSE) %>%
    ma$qc() %>%
    ma$normalize() %>%
    ma$annotate(summarize="hgnc_symbol")

# save result
save(expr, file=OUTFILE)
