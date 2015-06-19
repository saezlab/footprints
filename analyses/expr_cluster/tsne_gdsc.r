library(ggplot2)
gdsc = import('data/gdsc')
expr = gdsc$basal_expression()
tissues = unique(gdsc$tissues(minN=10))

pdf("tsne_gdsc.pdf", paper="a4r", width=26, height=20)
on.exit(dev.off)

for (tissue in tissues) {
    valid = intersect(colnames(expr), names(gdsc$tissues(tissue)))
    expr_sub = t(expr[,valid])
    result = bhtsneR::tsne(expr_sub, perplexity=nrow(expr_sub)/4)

    print(qplot(result[,1], result[,2]) + ggtitle(tissue))
}
