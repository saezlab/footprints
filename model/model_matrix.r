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
zscores = zobj$scores
index = zobj$index

# adjust object for linear modelling
inh = index$effect=="inhibiting"
zscores[,inh] = -zscores[,inh]
zscores = t(zscores)

# fit model to pathway perturbations
pathway = ar$mask(index$pathway) + 0
pathway[,"EGFR"] = pathway[,"EGFR"] + pathway[,"MAPK"] + pathway[,"PI3K"]
pathway[,"TNFa"] = pathway[,"TNFa"] + pathway[,"NFkB"]

index = c(as.list(index), list(zscores=zscores, pathway=pathway))
# ^^ combination between parent.frame()+expl.data <-fix
mod = st$lm(zscores ~ 0 + pathway, data=index, min_pts=30, atomic="pathway") %>%
    transmute(gene = zscores,
              pathway = sub("^pathway", "", term),
              zscore = estimate,
              p.value = p.value) %>%
    group_by(gene) %>%
    mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
    ungroup()

zfit = ar$construct(zscore ~ gene + pathway, data=mod)
pval = ar$construct(p.adj ~ gene + pathway, data=mod)

# filter zfit to only include top 100 genes per pathway
zfit[apply(pval, 2, z -> z > b$minN(z, 100))] = 0

# save resulting object
save(zfit, file=OUTFILE)
