b = import('base')
io = import('io')

# better: mapping yaml->zscore file
YAML = commandArgs(TRUE)[1] %or% "new_index/MAPK/E-GEOD-51212.yaml"
EXPR = commandArgs(TRUE)[2] %or% "normalized/E-GEOD-51212.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "zscores/MAPK/E-GEOD-51212.RData"

yaml = io$read_yaml(YAML, drop=FALSE)
expr = io$load(EXPR)

# for each contrast, compute z-scores
record2zscore = function(record, expr) {
    if (is.list(expr))
        expr = expr[[record$platform]]

    pd = Biobase::pData(expr)
    emat = Biobase::exprs(expr)
    mapping = setNames(rownames(pd), pd$Source.Name)
    sub_control = mapping[record$control]
    sub_perturbed = mapping[record$perturbed]

    # only use arrays that passed qc, fail if we discard whole record
    rec_both = c(record$control, record$perturbed)
    sub_both = c(sub_control, sub_perturbed)
    missing = is.na(sub_both)
    if (any(missing))
        warning("discarding ", paste(rec_both[missing], collapse=", "))
    sub_control = sub_control[!is.na(sub_control)]
    sub_perturbed = sub_perturbed[!is.na(sub_perturbed)]
    if (length(sub_control) < 2 || length(sub_perturbed) < 1)
        stop("record ", record$accession, " failed")

    # get expression values from source name
    control = emat[,sub_control, drop=FALSE]
    perturbed = emat[,sub_perturbed, drop=FALSE]

    # build models
    mean_control= apply(control, 1, mean)
    sd_control = apply(control, 1, sd)
    mean_perturbed = apply(perturbed, 1, mean)
    logFC = mean_perturbed - mean_control
    model = loess(sd_control ~ mean_control)
    zscores = logFC / predict(model, mean_perturbed)

    # create index column
    index = record
    index$control = paste(index$control, collapse=";")
    index$perturbed = paste(index$perturbed, collapse=";")

    list(zscores=zscores, dscores=logFC, index=index)
}

records = lapply(yaml, function(r) record2zscore(r, expr))
index = lapply(records, function(x) x$index) %>% do.call(rbind, .)
zscores = lapply(records, function(x) x$zscores) %>% do.call(cbind, .)
dscores = lapply(records, function(x) x$dscores) %>% do.call(cbind, .)
result = list(index=index, zscores=zscores, dscores=dscores)

# separate index file w/ metadata derived from yaml [preferred?]
save(result, file=OUTFILE)
