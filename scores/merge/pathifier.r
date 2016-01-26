library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
hpc = import('hpc')

#FIXME:
# it looks like pathifier is unable to calculated the score for individual
# samples (as SPIA or paradigm do).
# for this, we most likely have to move back to tissue-level model
# building which requires a rewrite of this routine again
# (status: deferred, want to get some GDSC-only approach working)

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/reactome.RData"
EXPR = commandArgs(TRUE)[2] %or% "../../util/expr_cluster/corrected_expr.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "pathifier.RData"

#' Calculates Pathifier scores for one sample vs tissue-matched normals
#'
#' @param sample    A character ID of the sample to compute scores for
#' @param expr      An expression matrix with [genes x samples]
#' @param index     A data.frame with at least the fields `id`, `tissue`
#'                  (TCGA tissue identifier), and `type` (normal vs tumor)
#' @param genesets  A list of character vectors corresponding to
#'                  gene sets (e.g. pathways for pathway scores)
sample2scores = function(sample, expr, index, genesets) {
    library(dplyr)
    io = import('io')
    pathifier = import_package('pathifier')

    sample_tissue = filter(index, id==sample)$tissue
    tissue_normals = index %>%
        filter(tissue == sample_tissue & grepl("[nN]ormal", type))

	data = cbind(expr[,sample,drop=FALSE], expr[,tissue_normals$id])

    result = pathifier$quantify_pathways_deregulation(
        data = data,
        allgenes = rownames(data),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = c(FALSE, rep(TRUE, nrow(tissue_normals))),
        attempts = 100,
        min_exp = -Inf,
        min_std = 0.4
    )
}

# load pathway gene sets and tissues
genesets = io$load(INFILE)
data = io$load(EXPR)
expr = data$expr
index = data$index

samples = filter(index, !grepl("[nN]ormal", type))$id

# run pathifier in jobs
result = hpc$Q(sample2scores, sample=samples,
               const=list(expr=expr, index=index, genesets=genesets),
               memory=8192, n_jobs=50) %>%
    setNames(samples) %>%
    ar$stack(along=1)

# save results
save(result, file=OUTFILE)
