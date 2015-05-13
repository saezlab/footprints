# point of this file:
# - use the zscores to create a linear model
library(dplyr)
b = import('base')
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
index = c(as.list(index), list(zscores=zscores)) # combination between parent.frame()+expl.data <-fix
mod = st$lm(zscores ~ 0 + pathway, data=index) %>%
    transmute(gene = zscores,
              pathway = sub("^pathway", "", term),
              zscore = estimate,
              p.value = p.value) %>%
    group_by(gene) %>%
    mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
    ungroup()

zfit = ar$construct(zscore ~ gene + pathway, data=mod)
pval = ar$construct(p.adj ~ gene + pathway, data=mod)

###
### old code begin
###
selectFeatures = function(X, min=NULL, max=NULL, n, discard.value=0) {
    if (!is.null(min) && !is.null(max))
        stop("only min OR max supported")

    select_mask = matrix(F, nrow=nrow(X), ncol=ncol(X), dimnames=dimnames(X))

    if (is.null(min)) # select by maximum
        keep = apply(max, 2, val -> val>=b$maxN(val, n))
    else # select by minimum
        keep = apply(min, 2, val -> val<=b$minN(val, n))

    X[!keep] = discard.value
    X
}
zfit = selectFeatures(zfit, min=pval, n=100)
###
### old code end
###

#zfit[pval>0.01] = 0

save(zfit, file="model_linear.RData")
