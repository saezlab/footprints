library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

# from: https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R
ks.test.2 <- function (x, y, ..., alternative = c("two.sided", "less", "greater"), exact = NULL, maxCombSize=10000) 
{
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 1L) 
        stop("not enough 'x' data")
    PVAL <- NULL
    if (is.numeric(y)) {
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        y <- y[!is.na(y)]
        n.x <- as.double(n)
        n.y <- length(y)
        if (n.y < 1L) 
            stop("not enough 'y' data")
        if (is.null(exact)) {
            exact <- (n.x * n.y < maxCombSize)
            if(!exact)
                warning(paste("P-value not computed exactly because",
                              "of combined sample size"))
        }
        METHOD <- "Two-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        n <- n.x * n.y/(n.x + n.y)
        w <- c(x, y)
        z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
        if (length(unique(w)) < (n.x + n.y)) {
            if (exact) {
                warning("cannot compute exact p-value with ties")
                exact <- FALSE
            }
            else warning("p-value will be approximate in the presence of ties")
            z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
            TIES <- TRUE
        }
        STATISTIC <- switch(alternative, two.sided = max(abs(z)), 
            greater = max(z), less = -min(z))

        edge <- which.max(abs(z))
        ES <- z[edge]
        
        nm_alternative <- switch(alternative, two.sided = "two-sided", 
            less = "the CDF of x lies below that of y", greater = "the CDF of x lies above that of y")
        if (exact && (alternative == "two.sided") && !TIES) 
            PVAL <- 1 - .Call(stats:::C_pSmirnov2x, STATISTIC, n.x, n.y)
    }
    else {
        if (is.character(y)) 
            y <- get(y, mode = "function", envir = parent.frame())
        if (!is.function(y)) 
            stop("'y' must be numeric or a function or a string naming a valid function")
        METHOD <- "One-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        if (length(unique(x)) < n) {
            warning("ties should not be present for the Kolmogorov-Smirnov test")
            TIES <- TRUE
        }
        if (is.null(exact)) 
            exact <- (n < 100) && !TIES
        x <- y(sort(x), ...) - (0:(n - 1))/n
        STATISTIC <- switch(alternative, two.sided = max(c(x, 
            1/n - x)), greater = max(1/n - x), less = max(x))
        if (exact) {
            PVAL <- 1 - if (alternative == "two.sided")
                result = tryCatch({
                .C(C_pkolmogorov2x, p = as.double(STATISTIC), 
                  as.integer(n), PACKAGE = "stats")$p
                }, warning = function(w) {
                    warning(w)
                }, error = function(e) {
                    .Call(C_pKolmogorov2x, STATISTIC, n)
                }, finally = {
                })

            else {
                pkolmogorov1x <- function(x, n) {
                  if (x <= 0) 
                    return(0)
                  if (x >= 1) 
                    return(1)
                  j <- seq.int(from = 0, to = floor(n * (1 - 
                    x)))
                  1 - x * sum(exp(lchoose(n, j) + (n - j) * log(1 - 
                    x - j/n) + (j - 1) * log(x + j/n)))
                }
                pkolmogorov1x(STATISTIC, n)
            }
        }
        nm_alternative <- switch(alternative, two.sided = "two-sided", 
            less = "the CDF of x lies below the null hypothesis", 
            greater = "the CDF of x lies above the null hypothesis")
    }
    names(STATISTIC) <- switch(alternative, two.sided = "D", 
        greater = "D^+", less = "D^-")
    if (is.null(PVAL)) {
        pkstwo <- function(x, tol = 1e-06) {
            if (is.numeric(x)) 
                x <- as.double(x)
            else stop("argument 'x' must be numeric")
            p <- rep(0, length(x))
            p[is.na(x)] <- NA
            IND <- which(!is.na(x) & (x > 0))
            if (length(IND))
                p[IND] <- tryCatch({
                    tryRes <- .C(stats:::C_pkstwo, length(x[IND]), p = x[IND], 
                             as.double(tol), PACKAGE = "stats")$p
                }, warning = function(w) {
                    warning(w)
                }, error = function(e) {
                    tryRes <- .Call(stats:::C_pKS2, p = x[IND], tol)
                }, finally = {
                })
            p
        }
        PVAL <- ifelse(alternative == "two.sided", 1 - pkstwo(sqrt(n) * 
            STATISTIC), exp(-2 * n * STATISTIC^2))
    }
    PVAL <- min(1, max(0, PVAL))
    RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, 
        method = METHOD, data.name = DNAME, ES = ES, edge = edge)
    class(RVAL) <- "htest"
    return(RVAL)
}

#' Function to calculate Kolmogorov-Smirnoff (GSEA) score for a single signature
#'
#' @param set     A character describing the gene set (element in sigs)
#' @param expr    The gene expression matrix [genes x samples]
#' @param sigs    The list of signatures
#' @return        Result for ks(expr[,sample], sigs[set])
ks = function(index, expr, sigs, ...) {
    message(index$id)
    ctl = expr[, index$control, drop=FALSE] %catch% return(NULL)
    ptb = expr[, index$perturbed, drop=FALSE] %catch% return(NULL)

    emat = cbind(ctl, ptb)
    colnames(emat) = c(rep("ctl", ncol(ctl)), rep("ptb", ncol(ptb)))
    design = narray::mask(colnames(emat)) + 0
    fit1 = limma::lmFit(emat, design)
    contrast = limma::makeContrasts("ptb-ctl", levels=design)
    fit2 = limma::contrasts.fit(fit1, contrast)
    fit3 = limma::eBayes(fit2)

    diff_expr = sort(fit3$t[,1])
    score = function(sig) {
        matched = sort(na.omit(match(sig, names(diff_expr))))
        re = ks.test.2(matched, (1:length(diff_expr))[-matched])
		re$ES
    }
    result = sapply(sigs, score)
}

if (is.null(module_name())) {
    INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/speed_matrix.RData"
    EXPR = commandArgs(TRUE)[2] %or% "../../data/expr.RData"
    OUTFILE = commandArgs(TRUE)[3] %or% "pathways_mapped/ks_speed_matrix.RData"

    # only filter when we didn't select manually
    if (grepl("mapped", OUTFILE)) {
        MIN_GENES = 1
        MAX_GENES = Inf
        job_size = 1
    } else {
        MIN_GENES = 5
        MAX_GENES = 500
        job_size = 20
    }

    # load required data
    speed = io$load(EXPR)
    index = speed$records
    expr = speed$expr
    genesets = io$load(INFILE)

    scores = mapply(ks, index=index, expr=expr,
                    MoreArgs = list(sigs=genesets)) %>%
        t()
#    scores[is.infinite(scores)] = 20 # all positive

    filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
    index = do.call(bind_rows, lapply(index, filter_index))

    save(scores, index, file=OUTFILE)
}
