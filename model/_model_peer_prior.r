library(peer)
library(modules)
library(dplyr)
io = import('io')
ar = import('array')
sg = import('sanger_robject')
lm = import('lm')
an = import('anova')
plt = import('plots')

# load speed dscores
dscores = io$data('SPEED-Data/SPEED2dmats')$rma_none

# set up prior
index = io$read_table("../SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
    select(id, pathway) %>%
    filter(id %in% colnames(dscores))

dscores = dscores[,index$id]
mod = lm$fit(t(dscores)~0+pathway, data=index) %>%
    lm$selectFeatures()

model=PEER() # initialize peer
PEER_setNk(model, 50)#ncol(dscores)) # number of hidden factors (just large enough)
PEER_setPhenoMean(model, dscores) # phenotype = expression data
#PEER_setCovariates(model, as.matrix(covariates)) # optional: set covariates
PEER_setTolerance(model, 1) # optional: set convergence parameters; default: 0.001
PEER_setVarTolerance(model, 0.1) # default: 0.00001
PEER_setSparsityPrior(model, mod$fit) # optional: apply prior
PEER_setPriorEps(model,1,10) # default: 0.1,10.
#PEER_setNmax_iterations(model, 10000) # optional: set number of steps; default: 1000
PEER_update(model) # calculate the model
#PEER_getResiduals(model) # diagnostics

# calc drug assocs w/ hidden factors
# optional/later: GO-validate them
#w = PEER_getW(model) # weights of different factors in different experiments
x = PEER_getX(model)
rownames(x) = rownames(dscores)

expr = sg$getBASAL_EXPRESSION()
ar$intersect(expr, x, along=1)
scores = t(expr) %*% x
scores = ar$map(scores, along=1, scale)

tissues = sg$getTissues()
Ys = sg$getDrugResponseForCellLines()
ar$intersect(scores, tissues, Ys, along=1)

assocs.pan = an$calcAssocs(Ys, scores, covariate=tissues, p.adjust="fdr")
print(plt$drawVolcano(assocs.pan, top=60, log='y', base.size=0.2))
