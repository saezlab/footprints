record2pathway = function(rec, exp, lookup) {
    spia = import('../../util/spia')
    b = import('base/operators')

    print(rec)

    rownames(exp) = lookup$entrezgene[match(rownames(exp), lookup$hgnc_symbol)]
    exp = limma::avereps(exp[!is.na(rownames(exp)),])

    result = spia$spia(exp[,rec$perturbed], exp[,rec$control], per_sample=TRUE,
                       pathids=spia$speed2kegg, verbose=TRUE)

    colnames(result) = spia$kegg2speed[colnames(result)]
    result

    # instead, if per_sample=FALSE
    # setNames(result, spia$kegg2speed[names(result)])
}

library(dplyr)
io = import('io')
ar = import('array')
hpc = import('hpc')

OUTFILE = "spia.RData"

# get index, expr data for test set
speed = io$load('../../data/expr.RData')
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
               memory=4096, n_jobs=5) #%>% setNames(names(records)) #%>% ar$stack(along=1)
#check if %>% ... required for per_sample=T/F

save(scores, index, file=OUTFILE)
