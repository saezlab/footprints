library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../model/model_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.RData"
TISSUES = import('../../config')$tcga$tissues

# load z vectors and TCGA expression
vecs = io$load(INFILE)$model
expr = tcga$rna_seq(TISSUES)

# calculate scores from expr and speed vectors
ar$intersect(vecs, expr, along=1)
scores = t(expr) %*% vecs %>%
    ar$map(along=1, scale)

# save scores
save(scores, file=OUTFILE)
