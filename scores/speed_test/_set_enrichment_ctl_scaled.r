library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "go.RData"

# load gene lists for pathways
genelist = io$load(INFILE)

# get index, expr data for test set
speed = io$load('../../data/expr.RData')
keep = sapply(speed$records, function(x) identical(x$exclusion, "test-set"))
index = speed$records[keep]
expr = speed$expr[keep]

doGSEA = function(index, expr, sigs) {
    gsea = import('../../util/gsea')
    re = gsea$runGSEA(expr=expr, sigs=sigs)
    mean_ctl = apply(re[index$control,], 2, mean)
    sd_ctl = apply(re[index$control,], 2, sd)
    colMeans((re[index$perturbed,] - mean_ctl) / sd_ctl)
}

# perform GSEA
scores = mapply(doGSEA, index=index, expr=expr,
                MoreArgs=list(sigs=genelist)) %>%
    t() %>%
    ar$map(along=1, scale)

#scores = hpc$Q(doGSEA, index=index, expr=expr, const = list(sigs=genelist),
#        memory = 1024, n_jobs = length(expr)) %>%
#    setNames(names(index)) %>%
#    ar$stack(along=1)

filter_index = function(x) x[! names(x) %in% c('control', 'perturbed', 'exclusion')]
index = lapply(index, filter_index) %>%
    bind_rows()

save(scores, index, file=OUTFILE)
