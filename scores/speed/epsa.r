library(dplyr)
library(magrittr)
b = import('base')
io = import('io')
st = import('stats')

#' Summarize each experiment series to mean control and mean perturbed
#'
#' @param recs   List of recods
#' @param exps   List of gene expression matrices
#' @param field  Field to subset ("control" or "perturbed")
#' @return       List with expression matrix per pathway
summarize_experiments = function(recs, exps, field) {
    mapply(function(r,e) {
        message(r$id)
        narray::map(e[,r[[field]]], along=2, mean) %catch% NULL
    }, r=recs, e=exps) %>%
        narray::stack(along=2)
}

#' Get ordered list of gene symbols for all EPSA pathways
#'
#' @param ctl   A matrix of mean control expression per experiment
#' @param ptb   A matrix of mean perturbed expression per experiment
#' @return      Named (HGNC symbols) character vector of median fold changes
epsa_pathway = function(ctl, ptb, index) {
    one_pathway = function(path) {
        message(path)
        cc = ctl[[path]]
        pp = ptb[[path]]
        pp = pp * narray::rrep(index$sign[match(colnames(pp), index$id)], nrow(pp))
        narray::intersect(cc, pp, along=1)
        paths = t(cbind(cc, pp))
        type = c(rep(0, ncol(cc)), rep(1, ncol(pp)))

        de_genes = st$lm(paths ~ 0 + type, hpc_args=list(n_jobs=5)) %>%
            arrange(p.value) %>%
            head(100) %>%
            arrange(-estimate)

        med_fc = pp - cc[,colnames(pp)]
        med_fc = narray::map(med_fc, along=2, function(x) median(x, na.rm=TRUE))
        ordered_genes = sort(med_fc, decreasing=TRUE)
        ordered_genes = ordered_genes[names(med_fc) %in% de_genes$paths]
    }

    split_pathway = function(x)
        narray::split(x, along=2, subsets=b$grep("^([^.]+)", colnames(x)))

    ctl = split_pathway(ctl)
    ptb = split_pathway(ptb)
    signatures = b$lnapply(names(ctl), one_pathway)
}

#' Score perturbation experiments according to EPSA
#'
#' @param sigs  Named (pathways) list of numeric vector (HGNC/median FC)
score_experiments = function(sigs, ctl, ptb) {
    mean_fc = ptb - ctl[rownames(ptb),colnames(ptb)]

    scores = sapply(sigs, function(s) cor(s, mean_fcs[names(s)])) %>%
        setNames(names(sigs))
}

if (is.null(module_name())) {
    EXPR = commandArgs(TRUE)[1] %or% "../../data/expr.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "epsa.RData"

    # load zscores, model building function, and expression for each experiment
    expr = io$load(EXPR)

    idx_remove = c("control", "perturbed")
    sign_lookup = setNames(c(1,-1), c("activating", "inhibiting"))
    index = lapply(expr$records, function(x) x[setdiff(names(x), idx_remove)]) %>%
        do.call(bind_rows, .) %>%
        mutate(sign = sapply(effect, function(x) sign_lookup[x]))

    ctl = summarize_experiments(expr$records, expr$expr, "control")
    ptb = summarize_experiments(expr$records, expr$expr, "perturbed")
    signatures = epsa_pathway(ctl, ptb, index)

    index = index[match(rownames(scores), index$id),]
    save(scores, index, file=OUTFILE)
}
