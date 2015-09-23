# point of this file:
# - use the zscores to create a linear model
library(dplyr)
io = import('io')
ar = import('array')
st = import('stats')

data = io$load('../data/expr.RData')
records = data$records
expr = data$expr
zscores = io$load('zscores.RData')

index = records %>%
    lapply(function(x) x[sapply(x, length) == 1]) %>%
    bind_rows() %>%
    filter(is.na(exclusion))

# subset expr and records to train set
records = records[index$id]
expr = expr[index$id]

# prepare for model building
zscores = zscores[,index$id]
inh = index$effect=="inhibiting"
zscores[,inh] = -zscores[,inh]
zscores = t(zscores)

# fit model to pathway perturbations
mod = st$lm(zscores ~ 0 + pathway, data=index, min_pts=100,
            hpc_args=list(n_jobs=10, memory=2048)) %>%
    transmute(gene = zscores,
              pathway = sub("^pathway", "", term),
              zscore = estimate,
              p.value = p.value) %>%
    group_by(gene) %>%
    mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
    ungroup()

zfit = ar$construct(zscore ~ gene + pathway, data=mod)
pval = ar$construct(p.adj ~ gene + pathway, data=mod)

# take record, expr and calculate pathway scores
cutoff_recall_test = function(zfit, pval, zcut, pcut, records, expr) {
    library(dplyr)
    ar = import('array')

    # optimise the cutoff parameters on the test set
    calc_scores = function(rec, exp, zfit, pval, zcut, pcut) {
        ar$intersect(zfit, pval, exp)

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
    concordance = mapply(con_fun, ar$split(scores, along=2, drop=TRUE), paths, reverse)
    sum(concordance) / length(concordance)

#    # do statistical test here
#    1.0 - phyper(q = sum(concordance) - 1,    # number of white balls drawn
#                 m = length(concordance)    , # numer of white in urn
#                 n = length(concordance) * (ncol(zfit)-1), # number of black in urn
#                 k = length(concordance))     # number of balls drawn
}

df = import('data_frame')
cut_df = df$create_index(zcut=seq(-2,2,0.1), pcut=c(1e-5, 1e-4, 1e-3, 1e-2, seq(0.05,0.5,0.05)),
                         args=list(zfit=zfit, pval=pval, records=records, expr=expr),
                         expand_grid=TRUE)
result = df$call(cut_df, cutoff_recall_test, hpc_args=list(n_jobs=50, memory=2048))

# save resulting object
save(zfit, file=OUTFILE)
