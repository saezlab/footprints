# for each contrast, compute z-scores
expr2zscore = function(rec, emat) {
    message(rec$id)

    # get expression values from source name
    control = emat[,rec$control, drop=FALSE]
    perturbed = emat[,rec$perturbed, drop=FALSE]

    # build speed models
    mean_control= apply(control, 1, mean)
    sd_control = apply(control, 1, sd)
    mean_perturbed = apply(perturbed, 1, mean)
    logFC = mean_perturbed - mean_control
    model = loess(sd_control ~ mean_control)
    
    logFC / predict(model, mean_perturbed)
}

library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% "expr.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "zscores.RData"

data = io$load(INFILE)
records = data$records
expr = data$expr

result = mapply(expr2zscore, rec=records, emat=expr, SIMPLIFY=FALSE)

zscores = lapply(result, function(x) x$zscores) %>%
    ar$stack(along=2)

# separate index file w/ metadata derived from yaml [preferred?]
save(zscores, file=OUTFILE)
