library(modules)
library(dplyr)
io = import('io')
ar = import('array')
an = import('anova')
lm = import('lm')
sg = import('sanger_robject')
plt = import('plots')

###
### load speed data, index
###
zscores = t(io$load('../SPEED-Data/SPEED2zmats.RData')$rma_none)
#dscores = t(io$load('../SPEED-Data/SPEED2dmats.RData')$rma_none)
index = io$read_table("../SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
    select(id, pathway, cells, GSE, GPL) %>%
    dplyr::filter(id %in% rownames(zscores))

zscores = zscores[index$id,]

###
### fit model to pathway perturbations
###
#zscores = dscores * abs(zscores)
mod = lm$fit(zscores~0+pathway, data=index)
zfit = mod$fit
pval = mod$pval
zfit = lm$selectFeatures(zfit, min=pval, n=100)

###
### load sanger data
###
Ys = sg$getDrugResponseForCellLines('IC50s') # or AUC
tissues = sg$getTissues(minN=5)
expr = t(sg$getBASAL_EXPRESSION())
ar$intersect(expr, tissues, Ys, along=1)

###
### calculate scores on panel
###
expr = expr[,intersect(rownames(zfit),colnames(expr))]
zfit = zfit[intersect(colnames(expr), rownames(zfit)),]
scores = expr %*% zfit %>%
    ar$map(along=1, base::scale) # pathway across experiments

###
### save pdf w/ pan-cancer & tissue specific
###
pdf("lm_control.pdf", paper="a4r", width=26, height=20)

assocs.pan = an$calcAssocs(Ys, scores, covariate=tissues, p.adjust="fdr")
print(plt$drawVolcano(assocs.pan, top=40, log='y', base.size=0.2))

Yf = sg$filterDrugResponse(Ys, tissues, top=0.1, abs=0, delta=2) # sensitive cell lines
assocs.tissue = an$calcAssocs(Yf, scores, subsets=tissues, p.adjust="fdr", stack=T)
print(plt$drawVolcano(assocs.tissue, top=40, p=0.2, log='y', base.size=2))

dev.off()
