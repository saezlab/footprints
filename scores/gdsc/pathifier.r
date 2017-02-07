library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')
hpc = import('hpc')
genesets = import('../../util/genesets')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/reactome.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathways_mapped/pathifier.RData"
MIN_GENES = 5
MAX_GENES = 500

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

# load pathway gene sets and tissues
expr = gdsc$basal_expression()
tissues = gdsc$tissues(minN=10)
expr = t(expr) # this should work with along=-1
ar$intersect(tissues, expr, along=1)
expr = t(expr)

genesets = io$load(INFILE) %>%
    genesets$filter_genesets(rownames(expr), MIN_GENES, MAX_GENES)

# make compatible to call with one set in above function
for (i in seq_along(genesets))
    genesets[[i]] = setNames(list(genesets[[i]]), names(genesets)[i])

sample_tissues = unique(tissues)

# run pathifier in jobs
result = hpc$Q(sample2scores, sample_tissue=sample_tissues, genesets=genesets,
               const=list(expr=expr, tissues=tissues),
               memory=8192, job_size=100, fail_on_error=FALSE, expand_grid=TRUE)

result[sapply(result, class) == "try-error"] = NULL #TODO: should work with NA
result = ar$stack(result, along=2)

# save results
save(result, file=OUTFILE)
