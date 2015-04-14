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
dscores = t(io$load('../SPEED-Data/SPEED2zmats.RData')$rma_none)
#dscores = t(io$load('../SPEED-Data/SPEED2dmats.RData')$rma_none)
index = io$read_table("../SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
    select(id, pathway, cells, GSE, GPL) %>%
    dplyr::filter(id %in% rownames(dscores))

dscores = dscores[index$id,]

###
### load sanger data
###
Ys = sg$getDrugResponseForCellLines('IC50s') # or AUC
tissues = sg$getTissues(minN=5)
expr = t(sg$getBASAL_EXPRESSION())
ar$intersect(expr, tissues, Ys, along=1)

###
### fit model to pathway perturbations
###
#dscores = dscores * abs(dscores)
coeffs = ar$mask(index$pathway)

#for (i in 1:50) {
#    mod = lm$fit(dscores~0+., data=coeffs)
    mod = lm$fit(dscores~0+pathway, data=index) #FIXME: why is this not the same??
    dfit = mod$fit
    pval = mod$pval
    colnames(dfit) = colnames(pval) = colnames(coeffs) #TODO: check this!
    dfit = lm$selectFeatures(dfit, min=pval, n=100)

    # cross talk on a sample-wise basis (this may be unstable)
    coeffs = dscores %*% dfit

    # cross talk on a pathway basis TODO:
    # ...

    # normalize coeffs? (-> 1 where perturbed ,fraction<=1 where else)
#}

###
### calculate scores on panel
###
expr = expr[,intersect(rownames(dfit),colnames(expr))]
dfit = dfit[intersect(colnames(expr), rownames(dfit)),]
scores = expr %*% dfit %>%
    ar$map(along=1, base::scale) # pathway across experiments

###
### save pdf w/ pan-cancer & tissue specific
###
pdf("lm_structure.pdf", paper="a4r", width=26, height=20)

assocs.pan = an$calcAssocs(Ys, scores, covariate=tissues, p.adjust="fdr")
print(plt$drawVolcano(assocs.pan, top=40, log='y', base.size=0.2))

Yf = sg$filterDrugResponse(Ys, tissues, top=0.1, abs=0, delta=2) # sensitive cell lines
assocs.tissue = an$calcAssocs(Yf, scores, subsets=tissues, p.adjust="fdr", stack=T)
print(plt$drawVolcano(assocs.tissue, top=40, p=0.2, log='y', base.size=2))

dev.off()
