io = import('io')
ar = import('array')
expr2scores = import('../../scores/speed/speed')$expr2scores

#' Model with single zscore, not consensus
#'
#' @param zdata  A list with the zscore matrix and index object
#' @return       The coefficients matrix [gene x experiment]
zdata2model = function(zdata) {
    b = import('base')

    index = zdata$index
    zscores = t(t(zdata$zscores) * index$sign)
    zscores[apply(zscores, 2, function(p) !b$min_mask(abs(p), 100))] = 0

    list(model = zscores)
}

if (is.null(module_name())) {
    # zs$zscores : z-scores genes x experiments
    # zs$index   : index df w/ id=experiment, other metadata
    zdata = io$load('../../data/zscores.RData')

    # expr$expr    : list[experiments] of expression matrices
    # expr$records : index as list
    expr = io$load('../../data/expr.RData')

    scores = clustermq::Q(expr2scores, id=zdata$index$id, job_size=50,
              const = list(expr=expr, zdata=zdata, zdata2model=zdata2model)) %>%
        setNames(zdata$index$id) %>%
        ar$stack(along=1) %>%
        ar$map(along=1, scale) # each sig should assign same total score

    index = dplyr::select(zdata$index, -exclusion)
    stopifnot(zdata$index$id == rownames(scores))

    # "pathways" are in cols
    save(scores, index, file="sigs_scores.RData")
}
