# point of this file:
# - use the different models model to calculate scores in
#   - input data
#   - ccle data
#   - tcga data (later)
b = import('base')
io = import('io')
ar = import('array')
ccle = import('data/ccle')

INFILE = commandArgs(TRUE)[1] %or% "../../model/model_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.RData"

# load zscores
zfit = io$load(INFILE)$model

# load sanger data
expr = t(ccle$basal_expression())
#ar$intersect(expr, tissues, along=1)

# calculate scores on panel
expr = expr[,intersect(rownames(zfit),colnames(expr))]
zfit = zfit[intersect(colnames(expr), rownames(zfit)),]
scores = expr %*% zfit %>%
    ar$map(along=1, base::scale) # pathway across experiments

save(scores, file=OUTFILE)
