b = import('base')
io = import('io')
ar = import('array')
spia = import('../../util/spia')
gdsc = import('data/gdsc')
hpc = import('hpc')

OUTFILE = commandArgs(TRUE)[1] %or% "spia.RData"

#' Calculates SPIA scores for one sample vs all other tissues
#'
#' @param sample   A character ID of the sample to compute scores for
#' @param expr     An expression matrix with [genes x samples]
#' @param tissues  A named (COSMIC ID) vector of TCGA tissues
#' @param spia     A loaded `spia` module to keep the paths form master
#'                 (this should not be required with zmq `hpc` module)
sample2scores = function(sample, expr, tissues, spia) {
    library(dplyr)

    sample_tissue = tissues[sample]
#   other_tissue = setdiff(names(tissues)[tissues == sample_tissue], sample)
    other_tissue = names(tissues)[tissues != sample_tissue]

    spia$spia(sample, other_tissue, data=expr, pathids=spia$speed2kegg)
}

# load pathway gene sets and tissues
expr = gdsc$basal_expression()
tissues = gdsc$tissues()
expr = t(expr) # this should work with along=-1
ar$intersect(tissues, expr, along=1)
expr = t(expr)

# map gene expression from HGNC to Entrez IDs
lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
    biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"),
    filter="hgnc_symbol", values=rownames(expr), mart=.)
rownames(expr) = lookup$entrezgene[match(rownames(expr), lookup$hgnc_symbol)]
expr = limma::avereps(expr[!is.na(rownames(expr)),])

# run spia in jobs and save
result = hpc$Q(sample2scores, sample=samples,
               const=list(expr=expr, tissues=tissues, spia=spia),
               memory=8192, n_jobs=30) %>%
    setNames(samples) %>%
    ar$stack(along=1) %>%
    ar$map(along=1, scale)

colnames(result) = spia$kegg2speed[colnames(result)]

save(result, file=OUTFILE)
