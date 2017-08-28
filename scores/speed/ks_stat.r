library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

#' Function to calculate Kolmogorov-Smirnoff (GSEA) score for a single signature
#'
#' @param set     A character describing the gene set (element in sigs)
#' @param expr    The gene expression matrix [genes x samples]
#' @param sigs    The list of signatures
#' @return        Result for ks(expr[,sample], sigs[set])
ks = function(index, expr, sigs, ...) {
    message(index$id)
    ctl = expr[, index$control]
    ptb = expr[, index$perturbed]

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
        re = ks.test(matched, (1:length(diff_expr))[-matched])
        unname(-log10(re$p.value) * sign(re$statistic))
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
        setNames(names(index)) %>%
        ar$stack(along=1)

    filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
    index = do.call(bind_rows, lapply(index, filter_index))

    save(scores, index, file=OUTFILE)
}
