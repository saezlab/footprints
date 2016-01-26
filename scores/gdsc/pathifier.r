import('base/operators')
io = import('io')
ar = import('array')
hpc = import('hpc')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/reactome.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathifier.RData"

pathifier_per_sample = function(sample, tissues, genesets, expr) {
    io = import('io')
    tcga = import('data/tcga')
    pathifier = import_package('pathifier')

    sample_tissue = tissues[sample]
    controls = setdiff(names(tissues)[tissues == sample_tissue], sample)

	data = cbind(expr[,sample,drop=FALSE], expr[,controls,drop=FALSE])
	is_normal = c(FALSE, rep(TRUE, length(controls)))

    result = pathifier$quantify_pathways_deregulation(
        data = data,
        allgenes = rownames(expr),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = is_normal,
#        attempts = 100, # default settings
        min_exp = -Inf, #TODO: std=4, applicable to voom?
#        min_std = 0.4
    )

    re = do.call(cbind, lapply(result$scores, c))
    rownames(re) = colnames(expr)
    re
}

# load pathway gene sets
genesets = io$load(INFILE)

tissues = gdsc$tissues(minN=10)
expr = t(gdsc$basal_expression())
ar$intersect(tissues, expr, along=1) #TODO: along=-1 should work
expr = t(expr)

samples = names(tissues) # COSMIC IDs

#TODO: remove this, update n_jobs
samples = sample(samples, 20)

# run pathifier in jobs
result = hpc$Q(pathifier_per_sample, sample=samples,
               const=list(tissues=tissues, genesets=genesets, expr=expr),
               memory=8192, n_jobs=10) #%>%
#    setNames(tissues) %>%
#    ar$stack(along=1)
#
## save results
#save(result, file=OUTFILE)
