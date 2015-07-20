library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../model/model_linear.RData"
EXPR = commandArgs(TRUE)[2] %or% "../../expr_cluster/corrected_expr.h5"
OUTFILE = commandArgs(TRUE)[3] %or% "speed.RData"

# load vectors
vecs = io$load(INFILE)
expr = t(io$h5load(EXPR, "/expr"))
ar$intersect(vecs, expr, along=1)

scores = t(expr) %*% vecs %>%
    ar$map(scale, along=1)

save(scores, file=OUTFILE)
