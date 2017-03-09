library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

#' Calculates the scores on each experiment excluding it in the signature
#'
#' @param id     A character ID of the current experiment to work on
#' @param expr   List of all input experiments, fields 'records' [list w/ id, etc.]
#'               and 'expr' [genes x arrays]
#' @param zdata  A list with fields 'index' providing experiment info
#'               [data.frame] and 'zscores' [genes x experiments]
#' @param zscore2model  The model building function (takes: zdata, hpc_args)
#' @return       Pathway scores for the current experiment
expr2scores = function(id, expr, zdata, zdata2model) {
    stopifnot(zdata$index$id == names(expr$records))
    ar = import('array')
    index = expr$records[[id]]
    expr = expr$expr[[id]]

    # build the model without the current experiment
    zdata$index = zdata$index[zdata$index$id!=id,]
    zdata$zscores = zdata$zscores[,colnames(zdata$zscores) != id]
    vecs = zdata2model(zdata)$model

    # calculate the scores for the current experiment
    ar$intersect(vecs, expr, along=1)
    mat = t(expr) %*% vecs
    ctl = mat[index$control,,drop=FALSE]
    ptb = mat[index$perturbed,,drop=FALSE]
    colMeans(ptb) - colMeans(ctl) # better w/o scale, but enough?
}

if (is.null(module_name())) {
    MODEL = commandArgs(TRUE)[1] %or% "../../model/model_matrix.r"
    EXPR = commandArgs(TRUE)[2] %or% "../../data/expr.RData"
    ZSCORES = commandArgs(TRUE)[3] %or% "../../data/zscores.RData"
    OUTFILE = commandArgs(TRUE)[4] %or% "speed_matrix.RData"

    # load zscores, model building function, and expression for each experiment
    zdata = io$load(ZSCORES)
    zdata2model = import_(sub("\\.r$", "", MODEL))$zscore2model
    expr = io$load(EXPR)

    scores = clustermq::Q(expr2scores, id=zdata$index$id, job_size=1,
              const = list(expr=expr, zdata=zdata, zdata2model=zdata2model)) %>%
        setNames(zdata$index$id) %>%
        ar$stack(along=1) %>%
        ar$map(along=1, scale) %>% # scale each pathway across all experiments
        ar$map(along=2, scale) # scale each experiment across all pathways

    index = dplyr::select(zdata$index, -exclusion)
    stopifnot(zdata$index$id == rownames(scores))

    save(scores, index, file=OUTFILE)
}
