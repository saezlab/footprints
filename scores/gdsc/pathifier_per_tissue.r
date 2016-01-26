library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/reactome.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathifier.RData"

#' Calculates Pathifier scores for one tissue vs all others
#'
#' @param sample    A character ID of the sample to compute scores for
#' @param expr      An expression matrix with [genes x samples]
#' @param index     A data.frame with at least the fields `id`, `tissue`
#'                  (TCGA tissue identifier), and `type` (normal vs tumor)
#' @param genesets  A list of character vectors corresponding to
#'                  gene sets (e.g. pathways for pathway scores)
sample2scores = function(sample, expr, tissues, genesets) {
    library(dplyr)
    io = import('io')
    pathifier = import_package('pathifier')

    sample_tissue = tissues[sample]
	other_tissue = setdiff(names(tissues)[tissues == sample_tissue] sample)

    result = pathifier$quantify_pathways_deregulation(
        data = expr,
        allgenes = rownames(expr),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = tissues != sample_tissue,
        attempts = 100,
        min_exp = -Inf,
        min_std = 0.01
    )
}

# load pathway gene sets and tissues
genesets = io$load(INFILE)
expr = gdsc$basal_expression()
tissues = gdsc$tissues()
expr = t(expr) # this should work with along=-1
ar$intersect(tissues, expr, along=1)
expr = t(expr)

# run pathifier in jobs
result = hpc$Q(sample2scores, sample=colnames(expr),
               const=list(expr=expr, tissues=tissues, genesets=genesets),
               memory=8192, n_jobs=length(sample_tissues)) %>%
    setNames(colnames(expr)) %>%
    ar$stack(along=1)

# save results
save(result, file=OUTFILE)
