library(modules)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
an = import('anova')
plt = import('plots')

if (is.null(module_name())) {
    factor_weight = io$load('models_4-8h.RData')[["50"]]$w

    index = io$read_table("../SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
        filter(id %in% rownames(factor_weight)) %>%
        select(id, pathway, cells)

    pathway = ar$mask(index$pathway[match(rownames(factor_weight), index$id)])
#    pathway = ar$mask(index$cells[match(rownames(factor_weight), index$id)])

    pdf("factor_assocs.pdf", paper="a4r", width=29.7, height=21)

    assoc.pan = an$calcAssocs(factor_weight, pathway, p.adjust="fdr")
    nassocs = sum(assoc.pan$pvalue < 0.05, na.rm=T)
    print(plt$drawVolcano(assoc.pan, top=60, log='y', base.size=2) +
          ggtitle(paste("covariate, factors =", 50, "assocs =", nassocs)))

#    assoc.tissue = an$calcAssocs(factor_weight, pathway, subsets=tissues, p.adjust="fdr")
#    nassocs = sum(assoc.tissue$pvalue < 0.05, na.rm=T)
#    print(plt$drawVolcano(assoc.tissue, top=60, log='y', base.size=2) +
#          ggtitle(paste("subsets, factors =", title, "assocs =", nassocs)))
#
    dev.off()
}
