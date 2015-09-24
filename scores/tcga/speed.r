library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../model/model_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.RData"

# load vectors
vecs = io$load(INFILE)

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - filter for last known alive? [where is this field?]
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

# calculate scores from expr and speed vectors
tt = tcga$tissues()
tissue2scores = function(t) {
    cur_vecs = vecs
    expr = tcga$rna_seq(t)
    ar$intersect(cur_vecs, expr, along=1)
    scores = t(expr) %*% cur_vecs %>%
        ar$map(along=1, scale)
}
scores = do.call(rbind, lapply(tt, tissue2scores))
scores = scores[!duplicated(rownames(scores)),]

# do lapply instead of for as soon as works, then ar$stack + save

save(scores, file=OUTFILE)
