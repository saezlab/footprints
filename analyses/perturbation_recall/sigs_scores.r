library(dplyr)
io = import('io')
ar = import('array')

#' Calculates the scores on each experiment excluding it in the signature
#'
#' This combines scores/speed/speed.r:expr2scores and model/model_matrix:zdata2scores
#'
#' @param id     A character ID of the current experiment to work on
#' @param expr   List of all input experiments, fields 'records' [list w/ id, etc.]
#'               and 'expr' [genes x arrays]
#' @param zdata  A list with fields 'index' providing experiment info
#'               [data.frame] and 'zscores' [genes x experiments]
#' @param zscore2model  The model building function (takes: zdata, hpc_args)
#' @return       Pathway scores for the current experiment
expr2scores = function(id, expr, zdata) {
    stopifnot(zdata$index$id == names(expr$records))
    b = import('base')
    ar = import('array')
    index = expr$records[[id]]
    expr = expr$expr[[id]]

    # build the model without the current experiment
    zdata$index = zdata$index[zdata$index$id!=id,]
    zdata$zscores = zdata$zscores[,colnames(zdata$zscores) != id]

    vecs = t(t(zdata$zscores) * zdata$index$sign)
    vecs[apply(vecs, 2, function(p) !b$min_mask(abs(p), 100))] = 0

    # calculate the scores for the current experiment
    ar$intersect(vecs, expr, along=1)
    mat = t(expr) %*% vecs
    ctl = mat[index$control,,drop=FALSE]
    ptb = mat[index$perturbed,,drop=FALSE]
    colMeans(ptb) - colMeans(ctl) # better w/o scale, but enough?
}

if (is.null(module_name())) {
    # zs$zscores : z-scores genes x experiments
    # zs$index   : index df w/ id=experiment, other metadata
    zdata = io$load('../../data/zscores.RData')

    # expr$expr    : list[experiments] of expression matrices
    # expr$records : index as list
    expr = io$load('../../data/expr.RData')

    scores = clustermq::Q(expr2scores, id=zdata$index$id, job_size=50,
              const = list(expr=expr, zdata=zdata)) %>%
        setNames(zdata$index$id) %>%
        ar$stack(along=1) %>%
        ar$map(along=1, scale) # do we want to sale this?

    index = dplyr::select(zdata$index, -exclusion)
    stopifnot(zdata$index$id == rownames(scores))

    # "pathways" are in cols
    save(scores, index, file="sigs_scores.RData")
}
