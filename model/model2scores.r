# point of this file:
# - use the different models model to calculate scores in
#   - input data
#   - gdsc data
#   - tcga data (later)
b = import('base')
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "model_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "scores_linear.RData"

# load zscores
zfit = io$load(INFILE)

# load sanger data
Ys = gdsc$getDrugResponse('IC50s') # or AUC
tissues = gdsc$getTissues(minN=5)
expr = t(gdsc$getBasalExpression())
ar$intersect(expr, tissues, Ys, along=1)

# calculate scores on panel
expr = expr[,intersect(rownames(zfit),colnames(expr))]
zfit = zfit[intersect(colnames(expr), rownames(zfit)),]
scores = expr %*% zfit #%>%
#    ar$map(along=1, base::scale) # pathway across experiments

save(scores, file=OUTFILE)
