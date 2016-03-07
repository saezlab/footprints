fit_model = function(index_orig, zscores_orig, resample=TRUE, myseed=NULL) {
    set.seed(myseed)

    library(dplyr)
    ar = import('array')
    st = import('stats')

    # resample zscores, index
    if (resample)
        index = index_orig %>%
            group_by(pathway) %>%
            sample_frac(1, replace=TRUE) %>%
            ungroup()
    else
        index = index_orig
    zscores = zscores_orig[index$id,]

############## this comes mostly from model_matrix script (refactor!)

    # fit model to pathway perturbations
    pathway = ar$mask(index$pathway) + 0
    pathway["EGFR",] = pathway["EGFR",] + pathway["MAPK",] + pathway["PI3K",]
    pathway["TNFa",] = pathway["TNFa",] + pathway["NFkB",]

    mod = st$lm(zscores ~ 0 + pathway, data=index, min_pts=30, atomic="pathway") %>%
        transmute(gene = zscores,
                  pathway = sub("^pathway", "", term),
                  zscore = estimate,
                  p.value = p.value) %>%
        group_by(gene) %>%
        mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
        ungroup()

    zfit = ar$construct(zscore ~ gene + pathway, data=mod)
    pval = ar$construct(p.adj ~ gene + pathway, data=mod)

    # filter zfit to only include top 100 genes per pathway
    zfit[apply(pval, 2, function(p) !b$min_mask(p, 100))] = 0

############## this comes from scores script (refactor!)
    io = import('io')
    EXPR = io$load('../../data/expr.RData')
    expr = EXPR$expr
    records = EXPR$records

    # calculate scores on full set
    expr2scores = function(index, expr, vecs) {
        ar$intersect(vecs, expr, along=1)
        mat = t(expr) %*% vecs
        ctl = mat[index$control,,drop=FALSE]
        ptb = mat[index$perturbed,,drop=FALSE]
        (colMeans(ptb) - colMeans(ctl)) #/ ar$map(ctl, along=1, sd)
    }

    scores = mapply(expr2scores, index=records, expr=expr,
        MoreArgs=list(vecs=zfit), SIMPLIFY=FALSE) %>%
        ar$stack(along=1) %>%
        ar$map(along=1, scale)
}

b = import('base', attach_operators=FALSE)
import('base/operators')
io = import('io')
hpc = import('hpc')

EXPR = commandArgs(TRUE)[1] %or% "../../data/expr.RData"
ZDATA = commandArgs(TRUE)[1] %or% "../../data/zscores.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "model_linear.RData"

# load speed data, index; filter for train set only
zdata = io$load(ZDATA)
index = zdata$index
zscores = t(zdata$zscores) * index$sign

re = fit_model(index_orig=index, zscores_orig=zscores, resample=TRUE, myseed=123)
#re = hpc$Q(fit_model, myseed=81925+1:100,
#           const=list(index_orig=index, zscores_orig=zscores, resample=TRUE),
#           memory=4096, n_jobs=100, log_worker=TRUE)

save(re, file="resample.RData")
