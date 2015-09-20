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

library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% "expr.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "scores.RData"

data = io$load(INFILE)
records = data$records
expr = data$expr

result = mapply(expr2zscore, rec=records, emat=expr, SIMPLIFY=FALSE)

index = lapply(result, function(x) x$index) %>%
    bind_rows()
### TEST SET TEST
# sample a third all accessions for each pathway, designate test set
set.seed(1829571)
test = index %>%
    select(pathway,accession) %>%
    distinct() %>%
    group_by(pathway) %>%
    do(sample_frac(.,0.3)) %>%
    ungroup()
excl = rep(NA, nrow(index))
excl[index$accession %in% test$accession] = "test-set"
index$exclusion = excl
### TEST SET TEST
zscores = lapply(result, function(x) x$zscores) %>%
    ar$stack(along=2)
dscores = lapply(result, function(x) x$dscores) %>%
    ar$stack(along=2) %>%
    impute::impute.knn()
dscores = dscores$data

# separate index file w/ metadata derived from yaml [preferred?]
save(index, zscores, dscores, file=OUTFILE)
