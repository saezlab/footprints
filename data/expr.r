#' Combine experiment index and normalized expression data
#'
#' This script will take as input all normalized gene expression data
#' and the experiment indices in .yaml files, and create an 'expr.RData'
#' object with two lists, one element for each accession-platform
#' combination:
#'
#' * expr: expression matrices of all samples that are used in at least
#'         one experiment
#' * records: content of the .yaml files; each element can have N
#'         experiments
#'
#' It will discard all experiments for which there are not at least
#' two controls and one perturbed arrays.
NULL

library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

filter_record = function(record, expr_avail) {
    record$control = intersect(record$control, expr_avail)
    record$perturbed = intersect(record$perturbed, expr_avail)
    if (length(record$control) > 1 && length(record$perturbed) > 0)
        record
    else
        NULL
}

subset_expr = function(id) {
    message(id)
    split_id = strsplit(id, "\\.")[[1]]

    if (!split_id[1] %in% names(yaml)) {
        warning("No index for expression object: ", id)
        return(NULL)
    }

    e = expr[[id]]
    yaml = yaml[names(yaml) == split_id[1]] # only one otherwise

    # convert sample identifiers to those used in AE/yaml
    pd = Biobase::pData(e)
    mapping = setNames(pd$Source.Name, rownames(pd))
    emat = Biobase::exprs(e)
    colnames(emat) = mapping[colnames(emat)]

    # filter out records that do not have enough samples available
    records = yaml[sapply(yaml, function(x) x$platform == split_id[2])] %>%
        lapply(function(r) filter_record(r, colnames(emat))) %>%
        b$omit$null()

    if (length(records) == 0) {
        warning("Discarding ID: ", id)
        return(NULL)
    }

    # subset expression matrix to the samples used
    in_both = sapply(records, function(x) c(x$control, x$perturbed)) %>%
        unlist() %>%
        unique()
    list(records=records, expr=emat[,in_both])
}

if (is.null(module_name())) {
    YAML = commandArgs(TRUE)[1] %or% list.files("new_index", "[0-9]+\\.yaml$",
                                     recursive=TRUE, full.names=TRUE)
    EXPR = commandArgs(TRUE)[2] %or% list.files("normalized", "\\.RData", full.names=TRUE)
    OUTFILE = commandArgs(TRUE)[3] %or% "./expr.RData"

    # load all index files
    yaml = lapply(YAML, function(y) io$read_yaml(y, drop=FALSE)) %>%
        do.call(c, .) %>%
        setNames(., lapply(., function(x) x$accession))

    # load all expression objects
    expr = lapply(EXPR, function(e) io$load(e)) %>%
        setNames(b$grep("/([^/]+)\\.RData", EXPR)) %>%
        do.call(c, .)

    # to check: discard expression objects that are not named properly
    keep = grepl("\\.", names(expr))
    if (any(!keep)) {
        warning("Discarding ", paste(names(expr)[!keep], collapse=", "))
        expr = expr[keep]
    }

    # get all data and corresponding index, save to file
    data = lapply(names(expr), subset_expr) %>%
        b$omit$null()
    records = lapply(data, function(x) x$records)
    expr = lapply(data, function(x) x$expr)
    save(records, expr, file=OUTFILE)
}
