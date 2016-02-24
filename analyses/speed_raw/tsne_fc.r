library(dplyr)
library(bhtsneR)
library(impute)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% "../../data/expr.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "tsne_fc.RData"

# calculate scores from expr and speed vectors
speed = io$load(INFILE)
index = speed$records
expr = speed$expr

# scaling: assume mean/sd across scores per sample is constant
# this protects against missing genes, etc in platform
expr2scores = function(index, expr) {
    mat = t(expr)
    ctl = mat[index$control,,drop=FALSE]
    ptb = mat[index$perturbed,,drop=FALSE]
    (colMeans(ptb) - colMeans(ctl)) #/ ar$map(ctl, along=1, sd)
}

scores = mapply(expr2scores, index=index, expr=expr, SIMPLIFY=FALSE) %>%
    ar$stack(along=1)

col_nas = colSums(!is.na(scores))
row_nas = rowSums(!is.na(scores))

scores = scores[row_nas/ncol(scores)>0.5,]
scores = scores[,col_nas/nrow(scores)>0.5]
scores = impute.knn(scores)$data

index = index[rownames(scores)]
index = lapply(index, function(x) x[!names(x) %in% c("control","perturbed","exclusion")]) %>%
	do.call(bind_rows, .)

dim2 = tsne(scores)
index$x = dim2[,1]
index$y = dim2[,2]

save(index, file=OUTFILE)
