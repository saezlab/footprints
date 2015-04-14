library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
lm = import('lm')
icgc = import('icgc')
an = import('anova')
plt = import('plots')

# construct pathway model from speed
dscores = io$load('../SPEED-Data/SPEED2dmats.RData')$rma_none

index = io$read_table("../SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
    select(id, pathway) %>%
    filter(id %in% colnames(dscores))

dscores = dscores[,index$id]
mod = lm$fit(t(dscores)~0+pathway, data=index)
vecs = lm$selectFeatures(mod$fit, min=mod$p, n=100)

# load clinical data, extract survival time+status
clinical = io$load('survival.RData')

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - filter for last known alive? [where is this field?]
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

# calculate scores from expr and speed vectors
expr = unique(clinical$icgc_sample_id) %>% icgc$getRNASeq(voom=TRUE)
ar$intersect(vecs, expr, along=1)
scores = t(expr) %*% vecs %>% ar$map(along=1, scale)

save(scores, file="expr_scores.RData")
