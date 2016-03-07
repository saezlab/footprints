library(dplyr)
library(reshape2)
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')

INFILE = "resample_matrix.RData"
OUTFILE = "score_matrix.RData"

zfits = io$load(INFILE)
expr = gdsc$basal_expression()

# calculate scores for all models and assemble to matrix
z2scores = function(z, expr) {
    ar$intersect(z, expr, along=1)
    re = t(expr) %*% z %>%
        ar$map(scale, along=1)
}
scores = lapply(zfits, function(x) z2scores(x, expr=expr)) %>%
    ar$stack(along=3)

save(scores, file=OUTFILE)
