b = import('base')
io = import('io')
ar = import('array')
spia = import('../../util/spia')
gdsc = import('data/gdsc')
hpc = import('hpc')

OUTFILE = commandArgs(TRUE)[1] %or% "spia.RData"
FILTER = as.logical(commandArgs(TRUE)[2]) %or% TRUE

#' Calculates SPIA scores for one sample vs all other tissues
#'
#' @param sample   A character ID of the sample to compute scores for
#' @param expr     An expression matrix with [genes x samples]
#' @param tissues  A named (COSMIC ID) vector of TCGA tissues
#' @param spia     A loaded `spia` module to keep the paths form master
#'                 (this should not be required with zmq `hpc` module)
#' @return         TODO
sample2scores = function(sample, expr, tissues, spia, pathids=NULL) {
#    b = import('base')
	spia = import('../../util/spia')
    sample_tissue = tissues[sample]
    other_tissue = setdiff(names(tissues)[tissues == sample_tissue], sample)
    spia$spia(sample, other_tissue, data=expr, pathids=pathids) #%catch% NA #TODO: shld b cvd by f_o_e
}

# load pathway gene sets and tissues
expr = gdsc$basal_expression()
tissues = gdsc$tissues(minN=10)
expr = t(expr) # this should work with along=-1
ar$intersect(tissues, expr, along=1)

expr = spia$map_entrez(t(expr))

if (FILTER) {
    genesets = spia$speed2kegg
} else {
    genesets = spia$pathids("hsa")
}

# make compatible to call with one set in above function
for (i in seq_along(genesets))
    genesets[[i]] = setNames(list(genesets[[i]]), names(genesets)[i])

# run spia in jobs and save
result = hpc$Q(sample2scores, sample=colnames(expr), pathids=genesets,
               const=list(expr=expr, tissues=tissues), expand_grid=TRUE,
               memory=8192, job_size=50, fail_on_error=FALSE)

result[sapply(result, class) == "try-error"] = NA
result = ar$stack(result, along=2)

if (FILTER)
    colnames(result) = spia$kegg2speed[colnames(result)]

save(result, file=OUTFILE)
