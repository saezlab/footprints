library(impute)
library(mclust)
io = import('io')

data = io$load('../../data/zscores.RData')
index = data$index
zscores = data$zscores

# number of genes per contrast
col_nas = colSums(!is.na(zscores))
hist(col_nas, 50)

# number of contrasts per number of genes
row_nas = rowSums(!is.na(zscores))
hist(row_nas, 50) # go enrichment about what we miss?

# subset where at least half exps have measurements, impute rest
zscores = zscores[row_nas/ncol(zscores)>0.5,]
zscores = zscores[,col_nas/nrow(zscores)>0.5]
zscores = impute.knn(zscores)$data

# see distribution of zscores themselves
save(zscores, index, file="zscores_filtered_imputed.RData")
