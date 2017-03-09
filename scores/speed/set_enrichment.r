library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

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

scores = clustermq::Q(doGSEA, index=index, expr=expr, const = list(sigs=genelist),
        memory = 1024, n_jobs = 50) %>%
    setNames(names(index)) %>%
    ar$stack(along=1) %>%
    ar$map(along=1, scale) #%>% # normalize different magnitude of pathways
#    ar$map(along=2, scale) # normalize total activation per experiment

filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
index = do.call(bind_rows, lapply(index, filter_index))

save(scores, index, file=OUTFILE)
