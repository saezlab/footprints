library(dplyr)
library(bhtsneR)
io = import('io')
ar = import('array')

#model = io$load('../../model/model_matrix.RData')
data = io$load('expr_filtered_imputed.RData')
expr = data$expr
index = data$index

#ar$intersect(expr, model, along=1)

scores = t(expr) #%*% model %>%
#    ar$map(along=1, scale) %>%
#   ar$map(along=2, scale)

stopifnot(rownames(scores) == index$id)
dim2 = tsne(scores)

index = index %>%
    mutate(x = dim2[,1],
           y = dim2[,2],
           pathway = ifelse(is.na(pathway), "basal", pathway),
           shape = as.factor(pathway == "basal")) %>%
    select(accession, platform, pathway, cells, treatment, hours, array, x, y, shape)

save(index, file="tsne_gex.RData")
