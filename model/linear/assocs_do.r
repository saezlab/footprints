library(modules)
b = import('base')
io = import('io')
ar = import('array')
an = import('anova')
lm = import('lm')
sg = import('sanger_robject')
plt = import('plots')
si = import('GDSC-Univariate/sangerInspect')
st = import('stats')
util = import('./assocs_util')

###
### load speed data, index
###
zscores = t(io$load('../SPEED-Data/SPEED2zmats.RData')$rma_none)
dscores = t(io$load('../SPEED-Data/SPEED2dmats.RData')$rma_none)
index = io$read_table("../SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
    select(pathway, cells, GSE, GPL) %>%
    dplyr::filter(id %in% rownames(zscores))

###
### fit model to pathway perturbations
###
#zscores = dscores * abs(zscores)
mod = lm$fit(zscores~0+pathway, data=index)
zfit = mod$fit
pval = mod$pval

###
### load sanger data
###
Ys = sg$getDrugResponseForCellLines('IC50s') # or AUC
tissues = sg$getTissues(minN=5)
expr = t(sg$getBASAL_EXPRESSION())
ar$intersect(expr, tissues, Ys, along=1)


###
### save pdf
###
pdf("do.pdf", paper="a4r", width=26, height=20)

###
### pan-cancer assocs
###
zfit = lm$selectFeatures(zfit, min=pval, n=100)
scores = util$calcScores(dscores, zfit)
util$drawHeatmap(scores, index)

#re = util$iterativeDiscard(zscores, dscores, index) # this kills most assocs
#util$drawHeatmap(re$scores, re$index, scale="column")
#scores2 = t(util$calcScores(expr, re$zfit)) # filtered zfit


scores2 = t(util$calcScores(expr, zfit, lm.fit=T)) # all arrays
save(scores2, file="pathways.RData")

assocs.pan = an$calcAssocs(Ys, scores2, covariate=tissues, p.adjust="fdr")
print(plt$drawVolcano(assocs.pan, top=40, log='y', base.size=0.2))

###
### tissue-specific drug assocs (do not use for now)
###
#Yf = sg$filterDrugResponse(Ys, tissues, top=0.1, abs=0, delta=2) # sensitive cell lines
##sig.drugs = apply(assocs.pan$pvalue, 1, function(x) any(x<0.05)) # significant pan assocs
##Yf = Yf[,sig.drugs]
#assocs.tissue = an$calcAssocs(Yf, scores, subsets=tissues, p.adjust="fdr", stack=T)
#print(plt$drawVolcano(assocs.tissue, top=40, p=0.2, log='y', base.size=2) + ggtitle(minN))

###
### plot fits for all significant pan-cancer assocs
###
#pindex = melt(assocs.pan$pvalue) %|>% x->filter(x,value<0.05) %|>% x->arrange(x, value)
pindex = melt(assocs.pan$pvalue) %>%
    dplyr::filter(value<0.05) %>%
    arrange(value)

for (i in 1:nrow(pindex))
    print(si$drugFit(Ys, pindex$Var1[i], scores2, pindex$Var2[i], tissues))

dev.off()
