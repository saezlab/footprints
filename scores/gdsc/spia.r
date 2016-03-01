b = import('base')
io = import('io')
ar = import('array')
spia = import('../../util/spia')
gdsc = import('data/gdsc')
hpc = import('hpc')

#' Calculates SPIA scores for one sample vs all other tissues
#'
#' @param sample   A character ID of the sample to compute scores for
#' @param expr     An expression matrix with [genes x samples]
#' @param tissues  A named (COSMIC ID) vector of TCGA tissues
#' @param spia     A loaded `spia` module to keep the paths form master
#'                 (this should not be required with zmq `hpc` module)
#' @return         TODO
sample2scores = function(sample, expr, tissues, spia, pathids=NULL) {
    sample_tissue = tissues[sample]
    other_tissue = setdiff(names(tissues)[tissues == sample_tissue], sample)
    spia$spia(sample, other_tissue, data=expr, pathids=pathids)
}

#' Maps row names from HGNC symbols to Entrez Gene IDs
#'
#' @param  expr  The expression matrix (genes x samples)
#' @return       An expression matrix with Entrez Gene IDs
map_entrez = function(expr) {
    # map gene expression from HGNC to Entrez IDs
    lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
        biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"),
        filter="hgnc_symbol", values=rownames(expr), mart=.)
    rownames(expr) = lookup$entrezgene[match(rownames(expr), lookup$hgnc_symbol)]
    expr = limma::avereps(expr[!is.na(rownames(expr)),])
}

if (is.null(module_name())) {
    OUTFILE = commandArgs(TRUE)[1] %or% "spia.RData"
    FILTER = as.logical(commandArgs(TRUE)[2]) %or% TRUE

    # load pathway gene sets and tissues
    expr = gdsc$basal_expression()
    tissues = gdsc$tissues(minN=10)
    expr = t(expr) # this should work with along=-1
    ar$intersect(tissues, expr, along=1)

    expr = map_entrez(t(expr))

    if (FILTER)
        pathids = spia$speed2kegg
    else
        pathids = NULL

    # run spia in jobs and save
    result = hpc$Q(sample2scores, sample=colnames(expr),
                   const=list(expr=expr, tissues=tissues, spia=spia, pathids=pathids),
                   memory=8192, job_size=20) %>%
        setNames(colnames(expr)) %>%
        ar$stack(along=1) %>%
        ar$map(along=1, scale)

    if (FILTER)
        colnames(result) = spia$kegg2speed[colnames(result)]

    save(result, file=OUTFILE)
}
