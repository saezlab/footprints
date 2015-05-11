# point of this file:
# - use the zscores to create a linear model
library(dplyr)
io = import('io')
ar = import('array')
st = import('stats')

# load speed data, index
zobj = io$load('../data/zscores.RData')
zscores = zobj$zscores
index = zobj$index

# adjust object for linear modelling
inh = index$effect=="inhibiting"
zscores[,inh] = -zscores[,inh]
zscores = na.omit(zscores) #TODO: handle this better
zscores = t(zscores)

# fit model to pathway perturbations
mod = st$lm(zscores ~ 0 + pathway, data=index) %>%
    transmute(gene = response,
              pathway = sub("^pathway", "", term),
              zscore = estimate,
              p.value = p.value) %>%
    group_by(gene) %>%
    mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
    ungroup()

zfit = ar$construct(zscore ~ gene + pathway, data=mod)
pval = ar$construct(p.adj ~ gene + pathway, data=mod)

#zfit = lm$selectFeatures(zfit, min=pval, n=100)
zfit[pval>0.01] = 0

save(zfit, file="model_linear.RData")
