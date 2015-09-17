record2pathway = function(rec, exp, lookup) {
    spia = import('../util/spia')
    b = import('base/operators')

    print(rec)

    rownames(exp) = lookup$entrezgene[match(rownames(exp), lookup$hgnc_symbol)]
    exp = limma::avereps(exp[!is.na(rownames(exp)),])

    result = spia$spia(exp[,rec$perturbed], exp[,rec$control], per_sample=FALSE,
                       pathids=spia$speed2kegg, verbose=TRUE)

    setNames(result, spia$kegg2speed[names(result)])
}

if (is.null(module_name())) {
    library(dplyr)
    io = import('io')
    ar = import('array')
    hpc = import('hpc')

    OUTFILE = "spia.RData"

    genesets = io$load("../util/genesets/reactome.RData")
    data = io$load('expr_subset.RData')
    expr = data$expr
    records = data$records

    # HGNC -> entrez gene lookup
    lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
        biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"), mart=.)

    scores = hpc$Q(record2pathway,
                   rec = records, exp = expr,
                   const = list(lookup = lookup),
                   memory=4096, n_jobs=5) %>% setNames(names(records)) %>% ar$stack(along=1)

    save(scores, file=OUTFILE)
}
