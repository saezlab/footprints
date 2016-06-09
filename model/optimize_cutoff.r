library(dplyr)
.ar = import('array')

calc_concordance = function(zfit, pval, zcut, pcut, records, expr) {
    # optimise the cutoff parameters on the test set
    calc_scores = function(rec, exp, zfit, pval, zcut, pcut) {
        .ar$intersect(zfit, pval, exp)

        # apply cutoffs
        zfit[pval > pcut] = 0
        if (zcut > 0)
            zfit[zfit < zcut] = 0
        else
            zfit[abs(zfit) < zcut] = 0

        # calculate mean scores, see if highest/lowest pathway matches record
        scores = t(zfit) %*% exp
        paths = rowMeans(scores[,rec$perturbed,drop=FALSE]) -
                         rowMeans(scores[,rec$control,drop=FALSE])

        list(pathway = rec$pathway,
             reverse = ifelse(rec$effect == "activating", FALSE, TRUE),
             scores = paths)
    }

    print(paste0("z: ", zcut, ", p: ", pcut))
    # get vector with TRUE if highest pathway is perturbed, FALSE otherwise
    # scale to factor out ||zfit * expr||
    re = mapply(calc_scores, rec=records, exp=expr,
        MoreArgs=list(zfit=zfit, pval=pval, zcut=zcut, pcut=pcut), SIMPLIFY=FALSE)
    paths = sapply(re, function(x) x$pathway)
    reverse = sapply(re, function(x) x$reverse)
    scores = sapply(re, function(x) x$scores) %>%
        ar$map(along=2, scale)

    # concordance: highest pathway score is the one activated (lowest/inhibited)
    con_fun = function(s,p,r) names(sort(s,decreasing=!r))[1] == p
    concordance = mapply(con_fun, .ar$split(scores, along=2, drop=TRUE), paths, reverse)
    sum(concordance) / length(concordance)

#    # do statistical test here
#    1.0 - phyper(q = sum(concordance) - 1,    # number of white balls drawn
#                 m = length(concordance)    , # numer of white in urn
#                 n = length(concordance) * (ncol(zfit)-1), # number of black in urn
#                 k = length(concordance))     # number of balls drawn
}

plot_concordance = function(concordance_result) {
    mat = .ar$construct(result ~ zcut + pcut, data=concordance_result)
    pheatmap::pheatmap(mat,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       display_numbers = TRUE)
}
