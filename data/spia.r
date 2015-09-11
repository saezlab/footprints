# todo: move somewhere else(?)

records2pathway = function(recs, exp, lookup) {
    spia = import('../util/spia')
    b = import('base/operators')

    record2pathway = function(rec, exp) {
        print(rec)

        rownames(exp) = lookup$entrezgene[match(rownames(exp), lookup$hgnc_symbol)]
        exp = limma::avereps(exp[!is.na(rownames(exp)),])

        spia$spia(exp[,rec$perturbed], exp[,rec$control], per_sample=FALSE,
                  pathids=spia$speed2kegg, verbose=TRUE) %catch% NA
    }

    lapply(recs, function(r) {
        result = record2pathway(r, exp)
        print(result)
        if (!identical(result, NA) && !identical(result, NULL)) # why NULL?
            setNames(result, spia$kegg2speed[names(result)])
        else
            NA
    })
}

if (is.null(module_name())) {
    library(dplyr)
    io = import('io')
    hpc = import('hpc')

    OUTFILE = "spia_scores.RData"

    genesets = io$load("../util/genesets/reactome.RData")
    data = io$load('expr.RData')
    expr = data$expr
    records = data$records

    # HGNC -> entrez gene lookup
    lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
        biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"), mart=.)

    scores = mapply(records2pathway, recs=records, exp=expr,
                    MoreArgs=list(lookup=lookup)) %>%
#    scores = hpc$Q(records2pathway,
#                   iter = list(recs = records, exp = expr),
#                   const = list(lookup = lookup),
#                   memory=4096, n_jobs=20) %>%
        unlist(recursive=FALSE) %>%
        unname()

    records = unname(unlist(records, recursive=FALSE))
    stopifnot(length(records) == length(scores))

    save(scores, records, file=OUTFILE)
}
