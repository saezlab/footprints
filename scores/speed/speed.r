library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% "../../model/model_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.RData"

# load vectors
vecs = io$load(INFILE)

# calculate scores from expr and speed vectors
speed = io$load('../../data/expr.RData')
keep = sapply(speed$records, function(x) identical(x$exclusion, "test-set"))
index = speed$records[keep]
expr = speed$expr[keep]

# scaling: assume mean/sd across scores per sample is constant
# this protects against missing genes, etc in platform
expr2scores = function(vecs, expr) {
    ar$intersect(vecs, expr, along=1)
    t(expr) %*% vecs %>% ar$map(along=1, scale)
}
scores = mapply(expr2scores, expr=expr, MoreArgs=list(vecs=vecs), SIMPLIFY=FALSE)

save(scores, index, file=OUTFILE)
