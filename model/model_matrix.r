# point of this file:
# - use the zscores to create a linear model
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')

INFILE = commandArgs(TRUE)[1] %or% "../data/zscores.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "model_matrix.RData"

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
pathway = ar$mask(index$pathway) + 0
pathway[,"EGFR"] = pathway[,"EGFR"] + pathway[,"MAPK"] + pathway[,"PI3K"]
#pathway[,"EGFR"] = 0.5*pathway[,"EGFR"] + 0.25*pathway[,"MAPK"] + 0.25*pathway[,"PI3K"]
mod = lm(zscores ~ 0 + pathway)
coeff = coef(summary(mod))
zfit = do.call(rbind, lapply(coeff, function(x) x[,'Estimate']))
pval = do.call(rbind, lapply(coeff, function(x) x[,'Pr(>|t|)']))
colnames(zfit) = colnames(pval) = sub("^pathway", "", colnames(zfit))
rownames(zfit) = rownames(pval) = sub("^Response ", "", rownames(zfit))

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

save(zfit, file=OUTFILE)
