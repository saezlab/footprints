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

zscores = mapply(expr2zscore, rec=records, emat=expr, SIMPLIFY=FALSE) %>%
    ar$stack(along=2)

idx_remove = c("control", "perturbed")
sign_lookup = setNames(c(1,-1), c("activating","inhibiting"))
index = lapply(records, function(x) x[setdiff(names(x), idx_remove)]) %>%
    bind_rows() %>%
    mutate(sign = sapply(effect, function(x) sign_lookup[x]))

stopifnot(colnames(zscores) == index$id)

# separate index file w/ metadata derived from yaml [preferred?]
save(zscores, index, file=OUTFILE)
