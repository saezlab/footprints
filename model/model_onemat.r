# point of this file:
# - use the zscores to create a linear model
library(dplyr)
b = import('base')
io = import('io')

#' Fits a linear model on Z-scores
#'
#' @param zdata  A list with the zscore matrix and index object
#' @return       The coefficients matrix [gene x pathway]
zscore2model = function(zdata, hpc_args=NULL) {
    b = import('base')
    ar = import('array')
    df = import('data_frame')
    st = import('stats')

    index = zdata$index
    zscores = t(zdata$zscores)

    # create indicator matrix which pathway is activated/inhibited
    pathway = index$sign * (ar$mask(index$pathway) + 0)
    pathway[,"MAPK"] = pathway[,"MAPK"] + pathway[,"EGFR"]
    pathway[,"NFkB"] = pathway[,"NFkB"] + pathway[,"TNFa"]
    # add EGFR>PI3K link here?

    # if using intercept (does this make a difference?)
#    pathway = cbind('(Intercept)'=1, pathway)

    mod = st$lm(zscores ~ 0 + pathway, min_pts=30, atomic="pathway",
                hpc_args=hpc_args) %>%
        transmute(gene = zscores,
                  pathway = sub("^pathway", "", term),
                  zscore = estimate,
                  p.value = p.value) %>%
#        filter(pathway != "(Intercept)") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))

    zfit = ar$construct(zscore ~ gene + pathway, data=mod, fill=0)
    pval = ar$construct(p.value ~ gene + pathway, data=mod, fill=1)

    # filter zfit to only include top 100 genes per pathway
    model = zfit
    model[apply(pval, 2, function(p) !b$min_mask(p, 100))] = 0

    list(assocs=mod, model=model)
}

if (is.null(module_name())) {
    ZDATA = commandArgs(TRUE)[1] %or% "../data/zscores.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "model_matrix.RData"

    # load speed data, index; filter for train set only
    zdata = io$load(ZDATA)
    result = zscore2model(zdata, hpc_args=list(n_jobs=20, memory=2048))

    # save resulting object
    save(result, file=OUTFILE)
}
