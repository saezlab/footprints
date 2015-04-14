library(modules)
library(pheatmap)
library(dplyr)
io = import('io')
ar = import('array')
sg = import('sanger_robject')
lm = import('lm')
an = import('anova')
plt = import('plots')

dscores = io$data('SPEED-Data/SPEED2dmats')$rma_none
index = io$read_table("../SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
    select(id, pathway) %>%
    filter(id %in% colnames(dscores)) %>%
    arrange(pathway)
#dscores = dscores[,index$id]

models = io$load('Nfactors_models.RData')
w = t(models[["100"]]$w[index$id,])
rownames(w) = 1:nrow(w)

rownames(index) = index$id # for pheatmaop matching
pheatmap(w, annotation=index['pathway'], cluster_cols=F, cluster_rows=F, scale="row")

#expr = sg$getBASAL_EXPRESSION()
#ar$intersect(expr, x, along=1)
#scores = t(expr) %*% x
#scores = ar$map(scores, along=1, scale)
#
#tissues = sg$getTissues()
#Ys = sg$getDrugResponseForCellLines()
#ar$intersect(scores, tissues, Ys, along=1)
