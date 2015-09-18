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
        attempts = 1000,
        min_exp = -Inf,
        min_std = 0.4
    )

    do.call(cbind, lapply(result$scores, c)) %>%
        ar$map(along=1, subsets=colnames(data), mean) %>%
        ar$map(along=1, function(x) x['perturbed']-x['control'])
}

library(dplyr)
io = import('io')
ar = import('array')
hpc = import('hpc')

OUTFILE = "pathifier.RData"

# load gene sets for pathways
genesets = io$load("../../util/genesets/reactome.RData")

# get index, expr data for test set
speed = io$load('../../data/expr.RData')
keep = sapply(speed$records, function(x) identical(x$exclusion, "test-set"))
index = speed$records[keep]
expr = speed$expr[keep]

scores = hpc$Q(record2pathway,
               rec = index, exp = expr,
               const = list(genesets=genesets),
               memory=1024, n_jobs=20) #%>%
#    setNames(names(records)) %>%
#    ar$stack(along=1)

save(scores, index, file=OUTFILE)
