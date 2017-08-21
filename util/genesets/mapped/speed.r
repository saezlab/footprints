library(magrittr)
b = import('base')
io = import('io')

INFILE = commandArgs(TRUE)[1] %or% "../../../model/model_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_matrix.RData"

sets = io$load(INFILE) %$%
    narray::collect(.$model != 0, along=1) %>%
    narray::split(along=2)

save(sets, file=OUTFILE)
