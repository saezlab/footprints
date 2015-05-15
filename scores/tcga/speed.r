library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
icgc = import('data/icgc')

INFILE = commandArgs(TRUE)[1] %or% "../../model/model_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed.RData"

# load vectors, clinical data
vecs = io$load(INFILE)
clinical = icgc$clinical()

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - filter for last known alive? [where is this field?]
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

# calculate scores from expr and speed vectors
expr = unique(clinical$icgc_sample_id) %>% icgc$rna_seq(voom=TRUE)
ar$intersect(vecs, expr, along=1)
scores = t(expr) %*% vecs %>% ar$map(along=1, scale)

save(scores, file=OUTFILE)
