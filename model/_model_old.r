# point of this file:
# - use the zscores to create a linear model
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
lm = import('./lm')

# load speed data, index
zscores = io$load('../../jsr-gdsc/SPEED-Data/SPEED2zmats.RData')$rma_none
index = io$read_table("../../jsr-gdsc/SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
    select(id, pathway, cells, GSE, GPL) %>%
    dplyr::filter(id %in% colnames(zscores))
zscores = t(zscores[,index$id])

## adjust object for linear modelling
#inh = index$effect=="inhibiting"
#zscores[,inh] = -zscores[,inh]
#zscores = na.omit(zscores) #TODO: handle this better
#zscores = t(zscores)

###
### begin old code
###
mod = lm$fit(zscores~0+pathway, data=index)
zfit = mod$fit
pval = mod$pval
###
### end old code
###

## fit model to pathway perturbations
#mod = st$lm(zscores ~ 0 + pathway, data=index) %>%
#    transmute(gene = response,
#              pathway = sub("^pathway", "", term),
#              zscore = estimate,
#              p.value = p.value) %>%
#    group_by(gene) %>%
#    mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
#    ungroup()
#
#zfit = ar$construct(zscore ~ gene + pathway, data=mod)
#pval = ar$construct(p.adj ~ gene + pathway, data=mod)

###
### old code begin
###
zfit = lm$selectFeatures(zfit, min=pval, n=100)
###
### old code end
###

#zfit[pval>0.01] = 0

save(zfit, file="model_old.RData")
