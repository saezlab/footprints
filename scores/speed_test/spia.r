record2pathway = function(rec, exp, lookup) {
    spia = import('../../util/spia')
    b = import('base/operators')

    print(rec)

    rownames(exp) = lookup$entrezgene[match(rownames(exp), lookup$hgnc_symbol)]
    exp = limma::avereps(exp[!is.na(rownames(exp)),])

    result = spia$spia(exp[,rec$perturbed], exp[,rec$control], per_sample=FALSE,
                       pathids=spia$speed2kegg, verbose=TRUE)

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
keep = sapply(speed$records, function(x) identical(x$exclusion, "test-set"))
index = speed$records[keep]
expr = speed$expr[keep]

# HGNC -> entrez gene lookup
lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
    biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"), mart=.)

#scores = mapply(record2pathway, rec=index, exp=expr,
#    MoreArgs=list(lookup=lookup), SIMPLIFY=FALSE)

scores = hpc$Q(record2pathway,
               rec = index, exp = expr,
               const = list(lookup = lookup),
               memory=4096, n_jobs=5) %>%
    setNames(names(index)) %>%
    ar$stack(along=1)

filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
index = lapply(index, filter_index) %>%
    bind_rows()

save(scores, index, file=OUTFILE)
