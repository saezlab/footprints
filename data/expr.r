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
invisible(NULL)

subset_expr = function(rec) {
    expr_id = paste(rec$accession, rec$platform, sep=".")
    e = expr[[expr_id]]

    if (is.null(e)) {
        warning("Discarding record (no expression available): ", rec$id)
        return(NULL)
    }

    # convert sample identifiers to those used in AE/yaml
    pd = Biobase::pData(e)
    mapping = setNames(pd$Source.Name, rownames(pd))
    emat = Biobase::exprs(e)
    colnames(emat) = mapping[colnames(emat)]

    # filter out records that do not have enough samples available
    rec$control = intersect(rec$control, colnames(emat))
    rec$perturbed = intersect(rec$perturbed, colnames(emat))
    if (length(rec$control) < 2 || length(rec$perturbed) < 1)
        warning("Discarding record (not enough samples): ", rec$id)
    else
        list(record=rec, expr=emat[,c(rec$control, rec$perturbed)])
}

if (is.null(module_name())) {
    library(dplyr)
    b = import('base')
    io = import('io')
    ar = import('array')

    EXPR = commandArgs(TRUE) %or% list.files("normalized", "\\.RData", full.names=TRUE)
    YAML = list.files("new_index", "[0-9]+\\.yaml$", recursive=TRUE, full.names=TRUE)
    OUTFILE = "expr.RData"

    # load all index files
    records = lapply(YAML, function(y) io$read_yaml(y, drop=FALSE)) %>%
        unlist(recursive=FALSE) %>%
        b$omit$null() %>%
        setNames(., lapply(., function(x) x$id))
    stopifnot(sum(duplicated(names(records))) == 0)

    # load all expression objects
    expr = lapply(EXPR, function(e) io$load(e)) %>%
        setNames(b$grep("/([^/]+)\\.RData", EXPR)) %>%
        unlist(recursive=FALSE)

    #TODO: check: discard expression objects that are not named properly
    keep = grepl("\\.", names(expr))
    if (any(!keep)) {
        warning("Discarding ", paste(names(expr)[!keep], collapse=", "))
        expr = expr[keep]
    }

    # get all data and corresponding index, save to file
    data = sapply(records, subset_expr, simplify=FALSE, USE.NAMES=TRUE) %>%
        b$omit$null()
    records = sapply(data, function(x) x$record, simplify=FALSE, USE.NAMES=TRUE)
### TEST SET TEST
# sample a third all accessions for each pathway, designate test set
set.seed(1829572)
test = records %>%
    lapply(function(x) x[c('pathway','accession')]) %>%
    bind_rows() %>%
    distinct() %>%
    group_by(pathway) %>%
    do(sample_frac(.,0.3)) %>%
    ungroup()
records = lapply(records, function(r) {
    if (r$accession %in% test$accession)
        r$exclusion = "test-set"
    else
        r$exclusion = NA
    r
})
message("length train set: ", sum(sapply(records, function(x) is.na(x$exclusion))))
message("length trest set: ", sum(sapply(records, function(x) identical(x$exclusion, "test-set"))))
### TEST SET TEST
    expr = sapply(data, function(x) x$expr, simplify=FALSE, USE.NAMES=TRUE)
    save(records, expr, file=OUTFILE)
}
