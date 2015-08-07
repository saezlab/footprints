library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
lincs = import('data/lincs')

INFILE = commandArgs(TRUE)[1] %or% "../../model/model_linear.RData"
INDEX = commandArgs(TRUE)[2] %or% "../../util/lincs_perturbation_qc/index.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "speed_linear.RData"

# load index, filter & load expression
index = unique(io$load(INDEX)$distil_id)
expr = lincs$get_z(cid=index, rid=lincs$projected, map.genes="hgnc_symbol")

# load vectors, calculate scores
vecs = io$load(INFILE)
ar$intersect(vecs, expr, along=1)

scores = t(expr) %*% vecs %>%
    ar$map(scale, along=1)

save(scores, file=OUTFILE)
