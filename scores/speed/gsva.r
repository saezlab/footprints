b = import('base')
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')
#gsea = import('../../util/gsea')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/go.RData"
EXPR = commandArgs(TRUE)[2] %or% "../../data/expr.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "pathways_mapped/gsea_go.RData"

# only filter when we didn't select manually
if (grepl("mapped", OUTFILE)) {
    MIN_GENES = 1
    MAX_GENES = Inf
    job_size = 1
} else {
    MIN_GENES = 5
    MAX_GENES = 500
    job_size = 20
}

# load required data
speed = io$load(EXPR)
index = speed$records
expr = speed$expr
genesets = io$load(INFILE)

#' Function to calculate GSVA score for a single signature
#'
#' @param set     A character describing the gene set (element in sigs)
#' @param expr    The gene expression matrix [genes x samples]
#' @param sigs    The list of signatures
#' @return        Result for GSVA(expr[,sample], sigs[set])
gsva = function(index, expr, sigs, ...) {
	mean_func = function(x) mean(x[index$perturbed]) - mean(x[index$control])
#	gsea$filter_genesets(rownames(expr), MIN_GENES, MAX_GENES)
    re = GSVA::gsva(expr=expr, gset.idx.list=sigs, parallel.sz=1, ...)$es.obs
    rowMeans(re[,index$perturbed,drop=FALSE]) - rowMeans(re[,index$control,drop=FALSE])
}

scores = hpc$Q(gsva, index=index, expr=expr, const = list(sigs=genesets),
        memory = 1024, n_jobs = 20) %>%
    setNames(names(index)) %>%
    ar$stack(along=1)

filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
index = do.call(bind_rows, lapply(index, filter_index))

save(scores, index, file=OUTFILE)
