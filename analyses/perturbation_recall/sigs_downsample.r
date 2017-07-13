#' Calculates the scores for all perturbation experiments
#'
#' @param expr   List of all input experiments, fields 'records' [list w/ id, etc.]
#'               and 'expr' [genes x arrays]
#' @param vecs   Pathway vectors matrix [genes x pathways]
#' @return       Pathway scores for the current experiment
expr2scores = function(expr, vecs) {
    # @param id     A character ID of the current experiment to work on
    id2scores = function(id, expr, vecs) {
        index = expr$records[[id]]
        expr = expr$expr[[id]]

        # calculate the scores for the current experiment
        narray::intersect(vecs, expr, along=1)
        mat = t(expr) %*% vecs
        ctl = mat[index$control,,drop=FALSE]
        ptb = mat[index$perturbed,,drop=FALSE]
        colMeans(ptb) - colMeans(ctl) # better w/o scale, but enough?
    }

    scores = lapply(names(expr), id2scores, expr=expr, vecs=vecs) %>%
        setNames(names(expr)) %>%
        narray::stack(along=2) %>%
        narray::map(along=1, scale)
}

#' Downsample the z-scores and recompute the model
#'
#' @param expr      List of expression matrice to score
#' @param zdata     A list with the zscore matrix and index object
#' @param frac      Fraction of experiments to sample (stratified by pathway)
#' @param hpc_args  Use HPC infrastructure
#' @return          The coefficients matrix [gene x experiment]
downsample_models = function(expr, zdata, frac=0.5, hpc_args=NULL) {
    zdata2model = import('../../model/model_matrix')$zscore2model
    library(dplyr)
    b = import('base')

    zdata$index = zdata$index %>%
        group_by(pathway) %>%
        sample_frac(frac) %>%
        ungroup()

    zdata$zscores = zdata$zscores[,zdata$index$id]
    model = zdata2model(zdata, hpc_args=hpc_args)$model
}

if (is.null(module_name())) {
    io = import('io')

    # zs$zscores : z-scores genes x experiments
    # zs$index   : index df w/ id=experiment, other metadata
    zdata = io$load('../../data/zscores.RData')

    # expr$expr    : list[experiments] of expression matrices
    # expr$records : index as list
    expr = io$load('../../data/expr.RData')

    frac = c(rep(0.5, 10), rep(0.25, 10))

#    models = clustermq::Q(downsample_models, frac=frac, job_size=1,
#              const = list(zdata = zdata,
#                           hpc_args = list(n_jobs=20)))]
    models = lapply(frac, downsample_models, zdata=zdata,
                    hpc_args = list(n_jobs=100))

#    scores = list(0.5 = lapply(scores[1:10], ),
#                  0.25 = scores[11:20])

    index = zdata$index

    # "pathways" are in cols
    save(scores, index, file="sigs_zscores.RData")
}
