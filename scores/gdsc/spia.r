b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
spia = import('../../util/spia')
gdsc = import('data/gdsc')

OUTFILE = commandArgs(TRUE)[1] %or% "pathways_mapped/spia.RData"
FILTER = as.logical(commandArgs(TRUE)[2]) %or% TRUE

#' Calculates SPIA scores for one sample vs all other tissues
#'
#' @param sample   A character ID of the sample to compute scores for
#' @param expr     An expression matrix with [genes x samples]
#' @param tissues  A named (COSMIC ID) vector of TCGA tissues
#' @return         TODO
sample2scores = function(sample, expr, tissues, pathids=NULL) {
	spia = import('../../util/spia')
    sample_tissue = tissues[sample]
    same_tissue = setdiff(names(tissues)[tissues == sample_tissue], sample)
    unname(spia$spia(sample, same_tissue, data=expr, pathids=pathids))
}

# load pathway gene sets and tissues
expr = gdsc$basal_expression()
tissues = gdsc$tissues(minN=10)
expr = t(expr) # this should work with along=-1
ar$intersect(tissues, expr, along=1)

expr = spia$map_entrez(t(expr))
if (FILTER) {
    pathids = spia$speed2kegg
} else {
    pathids = spia$pathids("hsa")
}

# run spia in jobs and save
result = b$expand_grid(sample = colnames(expr), pathids = pathids) %>%
    df$call(sample2scores, expr = expr, tissues = tissues,
            hpc_args = list(memory=4096, job_size=100, fail_on_error=FALSE))

result = ar$construct(result ~ sample + pathids, result)

if (FILTER)
    colnames(result) = spia$kegg2speed[colnames(result)]

save(result, file=OUTFILE)
