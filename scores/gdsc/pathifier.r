library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')
hpc = import('hpc')
gsea = import('../../util/gsea')

#' Calculates Pathifier scores for one tissue vs all others
#'
#' @param sample    A character ID of the sample to compute scores for
#' @param expr      An expression matrix with [genes x samples]
#' @param index     A data.frame with at least the fields `id`, `tissue`
#'                  (TCGA tissue identifier), and `type` (normal vs tumor)
#' @param genesets  A list of character vectors corresponding to
#'                  gene sets (e.g. pathways for pathway scores)
sample2scores = function(sample_tissue, expr, tissues, genesets) {
    library(dplyr)
    io = import('io')
    pathifier = import_package('pathifier')

    data = expr[, tissues==sample_tissue]

    # run for our pathways with default settings
    result = pathifier$quantify_pathways_deregulation(
        data = data,
        allgenes = rownames(expr),
        syms = genesets,
        pathwaynames = names(genesets),
        attempts = 100,
        min_exp = 4,
        min_std = 0.4
    )

	re = do.call(cbind, lapply(result$scores, c))
	rownames(re) = colnames(data)
	re
}

if (is.null(module_name())) {
    INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/reactome.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "pathifier.RData"
    MIN_GENES = 5
    MAX_GENES = 500

    # load pathway gene sets and tissues
    expr = gdsc$basal_expression()
    tissues = gdsc$tissues(minN=10)
    expr = t(expr) # this should work with along=-1
    ar$intersect(tissues, expr, along=1)
    expr = t(expr)

    genesets = io$load(INFILE) %>%
        gsea$filter_genesets(rownames(expr), MIN_GENES, MAX_GENES)

    sample_tissues = unique(tissues)

    # run pathifier in jobs
    result = hpc$Q(sample2scores, sample_tissue=sample_tissues,
                   const=list(expr=expr, tissues=tissues, genesets=genesets),
                   memory=8192, n_jobs=length(sample_tissues)) %>%
        ar$stack(along=1)

    # save results
    save(result, file=OUTFILE)
}
