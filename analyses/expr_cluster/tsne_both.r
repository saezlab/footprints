library(ggplot2)
io = import('io')

dset = io$load("../../gene_expr/corrected_expr.RData")
expr = dset$corrected
tissues = dset$covar

pdf("tsne_both.pdf", paper="a4r", width=26, height=20)
on.exit(dev.off)

for (tissue in unique(tissues)) {
    expr_sub = t(expr[,tissues==tissue])
    result = bhtsneR::tsne(expr_sub, perplexity=nrow(expr_sub)/4)

    print(qplot(result[,1], result[,2]) + ggtitle(tissue))
}
