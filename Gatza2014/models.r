#!/usr/bin/env Rscript
# Q: how much variance for each tissue-drug can be explained by tissue+<signature> features?

#' @param drug     drug index or identifer to subset \code{DRUGS}
#' @param tissue   tissue identifier to subset \code{DRUGS} and \code{FEATS}
#' @param DRUGS    drug matrix (cell lines x drugs)
#' @param FEATS    drug matrix (cell lines x features)
#' @param TISSUES  tissue vector (names = cell lines)
getMeanXValError = function(drug, tissue, DRUGS, FEATS, TISSUES) {
    stopifnot(rownames(DRUGS) == rownames(FEATS))
    stopifnot(rownames(DRUGS) == names(TISSUES))

    library(mlr)

    # subset data to current drug, tissue
    tissue_ = tissue==TISSUES
    tissue_[is.na(tissue_)] = FALSE
    feats_ = FEATS[tissue_,]
    target_ = DRUGS[tissue_, drug]

    # remove NA values
    discard = is.na(target_) | rowSums(is.na(feats_))>0
    feats_ = feats_[!discard,]
    target_ = target_[!discard]

    # make sure we have enough data points to build a model
    if (length(target_) < 11)
        return(NA)

    # construct data.frame that mlr can handle
    feats_ = as.data.frame(feats_)
    feats_$drug = target_
    colnames(feats_) = make.names(colnames(feats_))

    # train model, max 4 variables
    learner = makeLearner("regr.glmnet", dfmax=5) # nvars=4 instead?
    task = makeRegrTask(data=feats_, target='drug')
    result = crossval(learner, task, iters=10, models=TRUE)

    result$aggr
#    # report R^2 for the whole data set (not entirely right)
#    sapply(result$models, function(m)
#        summary(lm(truth ~ response, predict(m, task)))$r.squared
#    )
}

# drug response to 260 drugs (columns) <-> target
# pathway scores 19/50 pathways (columns) <-> features
# 19 different tissues (subset along rows)
if (is.null(module_name())) {
    io = import('io')
    ar = import('array')
    sg = import('sanger_robject')
    hpc = import('hpc')

    TISSUES = sg$getTissues()
    DRUGS = sg$getDrugResponseForCellLines()
    gatza = io$load('../GDSC-Scores/gatza.RData')[[1]]
    speed = io$load('../GDSC-Scores/speed.RData')[[1]]
    ar$intersect(TISSUES, DRUGS, gatza, speed)

    gatza = hpc$Q(getMeanXValError, #TODO: ar$stack results if value/vector/matrix?
        drug=colnames(DRUGS), tissue=unique(TISSUES),
        more.args=list(DRUGS=DRUGS, TISSUES=TISSUES, FEATS=gatza),
        expand.grid=TRUE, chunk.size=200, memory=512)

    speed = hpc$Q(getMeanXValError,
        drug=colnames(DRUGS), tissue=unique(TISSUES),
        more.args=list(DRUGS=DRUGS, TISSUES=TISSUES, FEATS=speed),
        expand.grid=TRUE, chunk.size=300, memory=512)

    index = expand.grid(drug=colnames(DRUGS), tissue=unique(TISSUES))
    index$gatza = unlist(gatza)
    index$speed = unlist(speed)
#    index$mean_gatza = sapply(gatza, mean)
#    index$sd_gatza = sapply(gatza, sd)
#    index$mean_speed = sapply(speed, mean)
#    index$sd_speed = sapply(speed, sd)

    save(index, file="models.RData") #TODO: load result, see how to process to save right thing
}
