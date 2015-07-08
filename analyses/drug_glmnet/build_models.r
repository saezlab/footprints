# -> just take combinations of the different scores
# -> train a model, report some cross-val error
# -> see which {,2-way} combination is the best {max 4 features}
# --
# -> or also: full model, how much variability can scores/mut/etc explain

.df = import('data_frame')
# Q: how much variance for each tissue-drug can be explained by tissue+<signature> features?

#' @param drug     drug index or identifer to subset \code{DRUGS}
#' @param tissue   tissue identifier to subset \code{DRUGS} and \code{FEATS}
#' @param DRUGS    drug matrix (cell lines x drugs)
#' @param FEATS    drug matrix (cell lines x features)
#' @param TISSUES  tissue vector (names = cell lines)
cv_err = function(formula, data=parent.frame(), min_pts=11, group=NULL, subsets=NULL, atomic=NULL, hpc_args=NULL) {
    #' @param ...   Arguments as defined in the data.frame row
    one_item = function(formula, data, subsets=NULL, ...) {
        args = list(...)
        # subset data according to subsets
        if (!is.null(subsets)) {
            data = lapply(data, function(x) x[subsets == args$subset,,drop=FALSE])
            args$subset = NULL
        }
        is_iterated = intersect(names(data), names(args))
        for (name in is_iterated)
            data[[name]] = data[[name]][, args[[name]], drop=TRUE]
        stopifnot(sapply(data[is_iterated], is.vector))

        # prepare ML model
        data = na.omit(as.data.frame(do.call(cbind, data)))
        colnames(data) = make.names(colnames(data)) # needed by mlr
        pts = nrow(data)
        if (pts < min_pts)
            return(NA)

        # train model, max 4 variables
        learner = mlr::makeLearner("regr.glmnet") #, dfmax=5) # nvars=4?
        task = mlr::makeRegrTask(data=data, target=as.character(formula[[2]]))
        result = mlr::crossval(learner, task, iters=10, models=FALSE, show.info=FALSE)
        result$aggr
    }

    idf = .df$from_formula(formula, data=data, group=group, subsets=subsets, atomic=atomic)
    .df$call(idf, one_item, hpc_args=hpc_args)
}


# drug response to 260 drugs (columns) <-> target
# pathway scores 19/50 pathways (columns) <-> features
# 19 different tissues (subset along rows)
if (is.null(module_name())) {
    io = import('io')
    ar = import('array')
    gdsc = import('data/gdsc')

    tissues = gdsc$tissues()
    drugs = gdsc$drug_response()
    gatza = io$load('../../scores/gdsc/gatza.RData')
    speed = io$load('../../scores/gdsc/speed_linear.RData')
    ar$intersect(tissues, drugs, gatza, speed, along=1)

    speedr = cv_err(drugs ~ speed, subsets=tissues, atomic="speed",
                    hpc_args = list(chunk.size=200, memory=512))
}
