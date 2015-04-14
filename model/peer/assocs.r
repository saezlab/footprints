library(modules)
io = import('io')
ar = import('array')
sg = import('sanger_robject')
an = import('anova')
plt = import('plots')

INFILE = commandArgs(TRUE)[1]
OUTFILE = commandArgs(TRUE)[2]

#TODO:? filter top 100 genes, recalc weights? (otherwise a lot of noise)

doAssocPlot = function(model, title, expr, tissues, Ys) {
    x = model$x
    ar$intersect(expr, x, along=1)
    scores = t(expr) %*% x
    scores = ar$map(scores, along=1, scale)
    ar$intersect(scores, tissues, Ys, along=1)

    assoc.pan = an$calcAssocs(Ys, scores, covariate=tissues, p.adjust="fdr")
    nassocs = sum(assoc.pan$pvalue < 0.05, na.rm=T)
    print(plt$drawVolcano(assoc.pan, top=60, log='y', base.size=0.2) +
          ggtitle(paste("covariate, factors =", title, "assocs =", nassocs)))

    assoc.tissue = an$calcAssocs(Ys, scores, subsets=tissues, p.adjust="fdr")
    nassocs = sum(assoc.tissue$pvalue < 0.05, na.rm=T)
    print(plt$drawVolcano(assoc.tissue, top=60, log='y', base.size=2) +
          ggtitle(paste("subsets, factors =", title, "assocs =", nassocs)))
}

models = io$load(INFILE)
expr = sg$getBASAL_EXPRESSION()
tissues = sg$getTissues()
Ys = sg$getDrugResponseForCellLines()

pdf(OUTFILE, paper="a4r", width=26, height=20)
for (i in seq_along(models))
    doAssocPlot(models[[i]], names(models)[i], expr, tissues, Ys)
dev.off()
