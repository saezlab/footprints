# point of this file:
# - use the different models model to calculate scores in
#   - input data
#   - gdsc data
#   - tcga data (later)

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
