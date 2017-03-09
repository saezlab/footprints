library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

#' Constructs a matrix of all control conditions
#'
#' @param data  A list with 'records' and 'expr' lists
#' @return      A matrix of genes x control arrays
all_control = function(data) {
    item2control = function(rec, emat) {
        control = emat[,rec$control, drop=FALSE]
        colnames(control) = sapply(colnames(control), function(x) paste(rec$id, x, sep="_"))
		control
    }
    mapply(item2control, rec=data$records, emat=data$expr, SIMPLIFY=FALSE) %>%
        ar$stack(along=2)
}

#' Compute z-scores for one contrast
#'
#' @param rec           The experiment record (`data.frame`, from the yaml files)
#' @param emat          The expression matrix [genes x experiments]
#' @param mean_control  The mean of the pooled control arrays
#' @param sd_control    The standard deviation of the pooled control arrays
expr2zscore = function(rec, emat, mean_control, sd_control) {
    message(rec$id)

    # get expression values from source name
    perturbed = emat[,rec$perturbed, drop=FALSE]

    # build speed models
    mean_perturbed = apply(perturbed, 1, mean)
    logFC = mean_perturbed - mean_control
    model = loess(sd_control ~ mean_control)

    logFC / predict(model, mean_perturbed)
}

#' Calculate Z-scores for all experiments
#'
#' @param data  A list with elements `records` and `expr`
#' @return      The `zscores` and `index` objects
data2zscores = function(data) {
    records = data$records
    expr = data$expr

	control = all_control(data)
    mean_control= apply(control, 1, function(x) mean(x, na.rm=TRUE))
    sd_control = apply(control, 1, function(x) sd(x, na.rm=TRUE))

    zscores = mapply(expr2zscore, rec=records, emat=expr,
            MoreArgs=list(mean_control=mean_control, sd_control=sd_control),
            SIMPLIFY=FALSE) %>%
        ar$stack(along=2)

    idx_remove = c("control", "perturbed")
    sign_lookup = setNames(c(1,-1), c("activating", "inhibiting"))
    index = lapply(records, function(x) x[setdiff(names(x), idx_remove)]) %>%
        do.call(bind_rows, .) %>%
        mutate(sign = sapply(effect, function(x) sign_lookup[x]))

    stopifnot(colnames(zscores) == index$id)

    list(zscores=zscores, index=index)
}

if (is.null(module_name())) {
    INFILE = commandArgs(TRUE)[1] %or% "expr.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "zscores_commonNULL.RData"

    data = io$load(INFILE)
	result = data2zscores(data)

    # separate index file w/ metadata derived from yaml [preferred?]
    zscores = result$zscores
    index = result$index
    save(zscores, index, file=OUTFILE)
}
