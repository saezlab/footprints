library(impute)
io = import('io')
ar = import('array')

data = io$load('expr.RData')
index = data$index
expr = data$expr

# number of genes per contrast
col_nas = colSums(!is.na(expr))
hist(col_nas, 50)
# number of contrasts per number of genes
row_nas = rowSums(!is.na(expr))
hist(row_nas, 50) # go enrichment about what we miss?

# subset where at least half exps have measurements, impute rest
expr = expr[row_nas/ncol(expr)>0.5,]
expr = expr[,col_nas/nrow(expr)>0.5]
expr = impute.knn(expr)$data

expr = t(expr)
ar$intersect(expr, index$id, along=1) # should probably support along=2 w/ df?
expr = t(expr)

# see distribution of expr themselves
save(expr, index, file="expr_filtered_imputed.RData")
