b = import('base')
io = import('io')
ar = import('array')
hpc = import('hpc')
gdsc = import('data/gdsc')

OUTFILE = commandArgs(TRUE)[1] %or% "spia.RData"

#' Calculate SPIA scores for each sample in each tissue
#'
#' @param sample   Current sample (character)
#' @param tissues  All tissues of samples in expr; needs same names as columns of expr
#' @param expr     Expression matrix [genes x samples]
spia_per_sample = function(sample, tissues, expr) {
    spia = import('../../util/spia')
    sample_tissue = tissues[sample]
    controls = setdiff(names(tissues)[tissues == sample_tissue], sample)
    spia$spia(sample, controls, data=expr, pathids=spia$speed2kegg)
}

tissues = gdsc$tissues(minN=10)
expr = t(gdsc$basal_expression())
ar$intersect(tissues, expr, along=1) #TODO: along=-1 should work
expr = t(expr)

samples = names(tissues) # COSMIC IDs

#FIXME: figure out why we don't get results for SPIA
#  might need to get TCGA normals, but I'd like to avoid that in order to be fair
#TODO: remove this, update n_jobs
samples = sample(samples, 20)

# run spia in jobs and save
result = hpc$Q(spia_per_sample, sample=samples,
               const=list(tissues=tissues, expr=expr),
               memory=8192, n_jobs=10) #%>%
#    ar$stack(along=1)

#colnames(result) = spia$kegg2speed[colnames(result)]

#save(result, file=OUTFILE)
