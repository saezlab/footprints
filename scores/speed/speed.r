library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% "../../model/model_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed.RData"

# load vectors
vecs = io$load(INFILE)

# calculate scores from expr and speed vectors
speed = io$load('../../data/dscores.RData')
expr = speed$scores
index = speed$index %>%
    select(-control, -perturbed)

ar$intersect(vecs, expr, along=1)
scores = t(expr) %*% vecs %>%
    ar$map(along=1, scale)

save(scores, index, file=OUTFILE)
