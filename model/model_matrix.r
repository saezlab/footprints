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

# filter zfit to only include top 100 genes per pathway
zfit[apply(pval, 2, z -> z > b$minN(z, 100))] = 0

# save resulting object
save(zfit, file=OUTFILE)
