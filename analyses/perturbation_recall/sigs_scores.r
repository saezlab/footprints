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
    library(dplyr)
    library(magrittr)
    stopifnot(zdata$index$id == names(expr$records))
    b = import('base')
    ar = import('array')
    index = expr$records[[id]]

    # subset control and perturbed matrices
    emat = expr$expr[[id]][,c(index$control, index$perturbed)]
    ctl = emat[,index$control]
    ptb = emat[,index$perturbed]

    # get top100 DE genes
    type = c(rep("ctl", ncol(ctl)), rep("ptb", ncol(ptb)))
    design = model.matrix(~ 0 + type)
    mod = limma::lmFit(emat, design)
    contrast = limma::makeContrasts("typeptb-typectl", levels=design)
    mod = limma::contrasts.fit(mod, contrast)
    mod = limma::eBayes(mod)
    top100 = as.data.frame(mod$p.value) %>%
        mutate(gene = rownames(.)) %>%
        arrange(`typeptb-typectl`) %>%
        head(100) %$%
        gene

    # get z vector for the top genes
    vec = zdata$zscores[,id,drop=FALSE][top100,,drop=FALSE]
    if (index$effect != "activating")
        vec = -vec

    # quantify signature in all other experiments
    others = expr$expr[,colnames(expr$expr) != id]
    ar$intersect(vec, emat, along=1)
    scores = t(emat) %*% vec
#    ctlmean = mean()
#    ptbmean = mean(t(ptb) %*% vec)
#    ptbmean - ctlmean
}

if (is.null(module_name())) {
    library(dplyr)
    io = import('io')
    ar = import('array')

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
