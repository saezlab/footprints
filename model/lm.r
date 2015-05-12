b = import('base')

# method: lm, rlm
fit = function(formula, data=NULL, method="lm") {
    if (is.matrix(data))
        data = as.data.frame(data)

    mod = lm(formula, data=data)
    coeff = coef(summary(mod))
    fit = do.call(rbind, lapply(coeff, function(x) x[,'Estimate']))
    pval = do.call(rbind, lapply(coeff, function(x) x[,'Pr(>|t|)']))

    varnames = all.vars(formula)
    if (length(varnames) == 2 && varnames[2] != '.')
        colnames(fit) = colnames(pval) = sub(paste0("^", varnames[2]), "", colnames(fit))

    rownames(fit) = rownames(pval) = sub("^Response ", "", rownames(fit))

    list(fit=fit, pval=pval)
}

# x: eg. expr
# vecs: eg. zfit
distance = function(x, vecs, lm.fit=F, scale=T) {
    expr = expr[,intersect(rownames(vecs),colnames(expr))]
    vecs = vecs[intersect(colnames(expr), rownames(vecs)),]

    if (lm.fit) { # assocs quite the same with both
        scores = sapply(1:nrow(expr), function(x) coef(lm(expr[x,]~0+vecs)))
        rownames(scores) = sub("vecs", "", rownames(scores))
        colnames(scores) = colnames(vecs)
    } else
        scores = t(expr %*% vecs)

    if (scale)
        scores = ar$map(scores, along=2, base::scale) # pathway across experiments
#        scores = ar$map(scores, along=1, base::scale) # experiment across pathways
        # along=2 does not change associations, but gives interpretable slopes
    scores
}

selectFeatures = function(X, min=NULL, max=NULL, n, discard.value=0) {
    if (!is.null(min) && !is.null(max))
        stop("only min OR max supported")

    select_mask = matrix(F, nrow=nrow(X), ncol=ncol(X), dimnames=dimnames(X))

    if (is.null(min)) # select by maximum
        keep = apply(max, 2, val -> val>=b$maxN(val, n))
    else # select by minimum
        keep = apply(min, 2, val -> val<=b$minN(val, n))

    X[!keep] = discard.value
    X
}
