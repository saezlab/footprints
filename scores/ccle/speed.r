b = import('base')
io = import('io')
ar = import('array')
ccle = import('data/ccle')

INFILE = commandArgs(TRUE)[1] %or% "../../model/model_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_matrix.RData"

# load model and expression data
zfit = io$load(INFILE)$model
expr = t(ccle$basal_expression())

# calculate scores on panel
expr = expr[,intersect(rownames(zfit),colnames(expr))]
zfit = zfit[intersect(colnames(expr), rownames(zfit)),]
scores = expr %*% zfit %>%
    ar$map(along=1, base::scale) # pathway across experiments

save(scores, file=OUTFILE)
