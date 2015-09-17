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
    zscores = logFC / predict(model, mean_perturbed)

    # create index column
    index = rec
    index$control = NULL
    index$perturbed = NULL

    list(zscores=zscores, dscores=logFC, index=index)
}

b = import('base')
io = import('io')
ar = import('array')

data = io$load('expr.RData')
records = data$records
expr = data$expr

result = mapply(expr2zscore, rec=records, emat=expr)

index = lapply(result, function(x) x$index) %>%
    bind_rows()
zscores = lapply(result, function(x) x$zscores) %>%
    ar$stack(along=1)
dscores = lapply(result, function(x) x$dscores) %>%
    ar$stack(along=1) %>%
    impute::impute.knn()
dscores = dscores$data

# separate index file w/ metadata derived from yaml [preferred?]
save(index, zscores, dscores, file=OUTFILE)
