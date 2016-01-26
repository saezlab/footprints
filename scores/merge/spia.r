b = import('base')
io = import('io')
ar = import('array')
spia = import('../../util/spia')
tcga = import('data/tcga')
gdsc = import('data/gdsc')
hpc = import('hpc')

EXPR = commandArgs(TRUE)[1] %or% "../../util/expr_cluster/corrected_expr.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "spia.RData"

#' Calculates SPIA scores for one sample vs tissue-matched normals
#'
#' @param sample  A character ID of the sample to compute scores for
#' @param expr    An expression matrix with [genes x samples]
#' @param index   A data.frame with at least the fields `id`, `tissue`
#'                (TCGA tissue identifier), and `type` (normal vs tumor)
#' @param spia    A loaded `spia` module to keep the paths form master
#'                (this should not be required with zmq `hpc` module)
sample2scores = function(sample, expr, index, spia) {
    library(dplyr)

    sample_tissue = filter(index, id==sample)$tissue
    tissue_normals = index %>%
        filter(tissue == sample_tissue & grepl("[nN]ormal", type))

    spia$spia(sample, tissue_normals$id, data=expr, pathids=spia$speed2kegg)
}

# load pathway gene sets and tissues
data = io$load(EXPR)
expr = data$expr
index = data$index
samples = filter(index, !grepl("[nN]ormal", type))$id

# map gene expression from HGNC to Entrez IDs
lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
    biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"),
    filter="hgnc_symbol", values=rownames(expr), mart=.)
rownames(expr) = lookup$entrezgene[match(rownames(expr), lookup$hgnc_symbol)]
expr = limma::avereps(expr[!is.na(rownames(expr)),])

# run spia in jobs and save
result = hpc$Q(sample2scores, sample=samples,
               const=list(expr=expr, index=index, spia=spia),
               memory=8192, n_jobs=100) %>%
    setNames(samples) %>%
    ar$stack(along=1) %>%
    ar$map(along=1, scale)

colnames(result) = spia$kegg2speed[colnames(result)]

save(result, file=OUTFILE)
