library(modules)
b = import('base')

calcScores = function(expr, zfit, lm.fit=F, scale=T) {
    expr = expr[,intersect(rownames(zfit),colnames(expr))]
    zfit = zfit[intersect(colnames(expr), rownames(zfit)),]

    if (lm.fit) { # assocs quite the same with both
        scores = sapply(1:nrow(expr), function(x) coef(lm(expr[x,]~0+zfit)))
        rownames(scores) = sub("zfit", "", rownames(scores))
        colnames(scores) = colnames(zfit)
    } else
        scores = t(expr %*% zfit)

    if (scale) {
        # along=2 does not change associations, but gives interpretable slopes
        # pathway across experiments FIXME?: f(#perturbed)
        scores = ar$map(scores, along=2, base::scale) # pathway across experiments
#        scores = ar$map(scores, along=1, base::scale) # experiment across pathways
    }
    scores
}

# lm.fit: use multiple regression?
# scale: arugment to pheatmap, "row" or "column"
drawHeatmap = function(scores, index, scale="none") {
    # re-index using order by pathway, then GPL, then GSE
    index = arrange(index, pathway, GPL, GSE)
    rownames(index) = index$subset # pheatmap needs this
    scores = scores[,index$subset]

    # plot heatmap with p-val filtered scores
    library(pheatmap)
    pheatmap(scores, cluster_cols=F, scale=scale,
             annotation=index[c('pathway','cells','GPL','GSE')])
}

iterativeDiscard = function(zscores, dscores, index) {
#    mod = apply(zscores, 2, z -> rlm(z~0+pathway, data=index, maxit=100))
#    coeff = lapply(mod, m -> coef(summary(m)))
#    zfit = do.call(rbind, lapply(coeff, function(x) x[,'Value']))
#    tval = do.call(rbind, lapply(coeff, function(x) x[,'t value']))
#    colnames(zfit) = colnames(tval) = sub("pathway", "", colnames(zfit))

    mod = lm(zscores~0+pathway, data=index)
    coeff = coef(summary(mod))
    zfit = do.call(rbind, lapply(coeff, function(x) x[,'Estimate']))
    pval = do.call(rbind, lapply(coeff, function(x) x[,'Pr(>|t|)']))
    colnames(zfit) = sub("pathway", "", colnames(zfit))
    rownames(zfit) = sub("Response ", "", rownames(zfit))

    scores = dscores %*% zfit
    bscores = ar$map(scores, along=2, x -> x==max(x))

    bpaths = mat$factorToModelMatrix(as.factor(index$pathway))
    keep = as.logical(rowSums(bscores & bpaths))
    print(sum(keep))

    if (all(keep))
        list(zfit=zfit, scores=t(scores), index=index, pval=pval)
    else
        iterativeDiscard(zscores[keep,], dscores[keep,], index[keep,])
}

# compile p-values for different tissues here, plot
Fisher.test = function(p) {
  Xsq = -2*sum(log(p))
  p.val = pchisq(Xsq, df = 2*length(p),lower.tail=FALSE)
  return(c(Xsq = Xsq, p.value = p.val))
}

