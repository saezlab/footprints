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
#model = "gene ~ EGFR + MAPK + PI3K + p53 + H2O2 + Hypoxia"
#fit = lavaan::sem(model, data = data)


gene = zscores[,2]
dfac = index$pathway
dmat = model.matrix(~0+index$pathway)
colnames(dmat) = sub("index\\$pathway", "", colnames(dmat))

#summary(lm(gene ~ 0+dmat))
summary(lm(gene ~ 0+dfac)) # same result

#summary(lavaan::sem("gene ~ dfac")) # 0 observations, same for dmat
#summary(lavaan::sem("gene ~ EGFR + H2O2 + Hypoxia + MAPK + p53 + PI3K + Trail", data=cbind(gene=gene, dmat)))

gmat = zscores[,sample(1:ncol(zscores),100)]
m = sapply(colnames(gmat), function(n) paste(n, "~ EGFR + H2O2 + Hypoxia + MAPK + p53 + PI3K + Trail"))
mod =  paste(c("MAPK ~ EGFR", m, sep="\n"))
summary(lavaan::sem(mod, data=cbind(gmat, dmat)))#, fit.measures=TRUE)

#TODO: run this for all genes, make plot + compare different options

# will likely have to use openmx because lavaan too slow
# -> try to get same example w/ 100 genes running + compare times
