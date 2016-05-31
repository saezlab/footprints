library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/go.RData"
EXPR = commandArgs(TRUE)[2] %or% "../../data/expr.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "go.RData"

# load required data
speed = io$load(EXPR)
index = speed$records
expr = speed$expr
genelist = io$load(INFILE)

doGSEA = function(index, expr, sigs) {
    library(dplyr)
    ar = import('array')
    gsea = import('../../util/gsea')
    mean_func = function(x) mean(x[index$perturbed]) - mean(x[index$control])

    gsea$runGSEA(expr=expr, sigs=sigs, normalize=TRUE) %>%
        ar$map(along=1, mean_func)
}

scores = hpc$Q(doGSEA, index=index, expr=expr, const = list(sigs=genelist),
        memory = 1024, n_jobs = 50) %>%
    setNames(names(index)) %>%
    ar$stack(along=1)

filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
index = lapply(index, filter_index) %>%
    bind_rows()

save(scores, index, file=OUTFILE)
