io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
plt = import('plot')

# load speed data, index
expr = io$load('../data/zscores.RData')
zscores = expr$zscores
index = expr$index
#TODO: inhibiting zscores = -1*zscores (or leave them out for now)

# fit model to pathway perturbations
mod = st$lm(zscores~0+pathway, data=index)
zfit = lm$selectFeatures(zfit, min=pval, n=100)

# load sanger data
Ys = gdsc$getDrugResponse('IC50s') # or AUC
tissues = sg$getTissues(minN=5)
expr = t(sg$getBasalExpression())
ar$intersect(expr, tissues, Ys, along=1)

# calculate scores on panel
expr = expr[,intersect(rownames(zfit),colnames(expr))]
zfit = zfit[intersect(colnames(expr), rownames(zfit)),]
scores = expr %*% zfit %>%
    ar$map(along=1, base::scale) # pathway across experiments

# save pdf w/ pan-cancer & tissue specific
pdf("lm_control.pdf", paper="a4r", width=26, height=20)

assocs.pan = an$calcAssocs(Ys, scores, covariate=tissues, p.adjust="fdr")
print(plt$drawVolcano(assocs.pan, top=40, log='y', base.size=0.2))

Yf = sg$filterDrugResponse(Ys, tissues, top=0.1, abs=0, delta=2) # sensitive cell lines
assocs.tissue = an$calcAssocs(Yf, scores, subsets=tissues, p.adjust="fdr", stack=T)
print(plt$drawVolcano(assocs.tissue, top=40, p=0.2, log='y', base.size=2))

dev.off()
