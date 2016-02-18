record2pathway = function(rec, exp, lookup) {
    spia = import('../../util/spia')
    b = import('base/operators')

    print(rec)

    rownames(exp) = lookup$entrezgene[match(rownames(exp), lookup$hgnc_symbol)]
    exp = limma::avereps(exp[!is.na(rownames(exp)),])

    result = spia$spia(exp[,rec$perturbed], exp[,rec$control], pathids=spia$speed2kegg)

#    # if per_sample=TRUE
#    colnames(result) = spia$kegg2speed[colnames(result)]
#    result

    setNames(result, spia$kegg2speed[names(result)])
}

library(dplyr)
import('base/operators')
io = import('io')
ar = import('array')
hpc = import('hpc')

EXPR = commandArgs(TRUE)[1] %or% '../../data/expr.RData'
OUTFILE = commandArgs(TRUE)[2] %or% "spia.RData"

# get index, expr data for test set
speed = io$load(EXPR)

# HGNC -> entrez gene lookup
lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
    biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"), mart=.)

#scores = mapply(record2pathway, rec=index, exp=expr,
#    MoreArgs=list(lookup=lookup), SIMPLIFY=FALSE)

result = hpc$Q(record2pathway,
               rec = speed$records, exp = speed$expr,
               const = list(lookup = lookup),
               memory=4096, n_jobs=10, fail_on_error=FALSE) %>%
    setNames(names(index))

errors = sapply(result, function(r) class(r) == "try-error")
if (any(errors)) {
    print(result[errors])
    result[errors] = NA
}
scores = ar$stack(result, along=1)

filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
index = lapply(speed$records[rownames(scores)], filter_index) %>%
    do.call(bind_rows, .)

save(scores, index, file=OUTFILE)
