library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/go.RData"
EXPR = commandArgs(TRUE)[2] %or% "../../data/expr.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "go.RData"

# load gene lists for pathways
genelist = io$load(INFILE)

# get index, expr data for test set
speed = io$load(EXPR)
keep = sapply(speed$records, function(x) identical(x$exclusion, "test-set"))
index = speed$records[keep]
expr = speed$expr[keep]

doGSEA = function(index, expr, sigs) {
    library(dplyr)
    ar = import('array')
    gsea = import('../../util/gsea')

    gsea$runGSEA(expr=expr, sigs=sigs, normalize=TRUE) %>%
        ar$map(along=1, function(x) mean(x[index$perturbed]) - mean(x[index$control]))
}

# perform GSEA
#scores = mapply(gsea$runGSEA, expr=expr,
#    MoreArgs=list(sigs=genelist, normalize=TRUE))

scores = hpc$Q(doGSEA, index=index, expr=expr, const = list(sigs=genelist),
        memory = 1024, n_jobs = length(expr)) %>%
    setNames(names(index)) %>%
    ar$stack(along=1)

filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
index = lapply(index, filter_index) %>%
    bind_rows()

save(scores, index, file=OUTFILE)
