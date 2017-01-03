record2pathway = function(rec, exp, genesets) {
    library(dplyr)
    import('base/operators')
    ar = import('array')
    pathifier = import_package('pathifier')

    print(rec)

    data = exp[,c(rec$control, rec$perturbed)]
    colnames(data) = c(rep("control", length(rec$control)),
                       rep("perturbed", length(rec$perturbed)))

    result = pathifier$quantify_pathways_deregulation(
        data = data,
        allgenes = rownames(data),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = colnames(data) == "control",
# default values
#        attempts = 100, # maybe set this higher and see if fewer NAs
#        min_exp = 4,
#        min_std = 0.4
    )

    do.call(cbind, lapply(result$scores, c)) %>%
        ar$map(along=1, subsets=colnames(data), mean) %>%
        ar$map(along=1, function(x) x['perturbed']-x['control'])
}

library(dplyr)
import('base/operators')
io = import('io')
ar = import('array')
hpc = import('hpc')

GENESETS = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/reactome.RData"
EXPR = commandArgs(TRUE)[2] %or% '../../data/expr.RData'
OUTFILE = commandArgs(TRUE)[3] %or% "pathifier.RData"

# load gene sets for pathways
genesets = io$load(GENESETS)

# get index, expr data for test set
speed = io$load(EXPR)
index = speed$records
expr = speed$expr

result = hpc$Q(record2pathway,
               rec = index, exp = expr,
               const = list(genesets=genesets),
               memory=1024, n_jobs=50, fail_on_error=FALSE) %>%
    setNames(names(index))

errors = sapply(result, function(r) class(r) == "try-error")
if (any(errors)) {
    print(result[errors])
    result[errors] = NA
}
scores = ar$stack(result, along=1)

filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
index = do.call(bind_rows, lapply(index, filter_index))

save(scores, index, file=OUTFILE)
