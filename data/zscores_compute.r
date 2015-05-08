b = import('base')
io = import('io')

# better: mapping yaml->zscore file
YAML = commandArgs(TRUE)[1] %or% "new_index/MAPK/E-GEOD-51212.yaml"
EXPR = commandArgs(TRUE)[2] %or% "normalized/E-GEOD-51212.RData"
ZFILE = commandArgs(TRUE)[3] %or% "zscores/MAPK/E-GEOD-51212.RData"

yaml = io$read_yaml(YAML, drop=FALSE)
expr = io$load(EXPR)

# for each contrast, compute z-scores
record2zscore = function(record, expr) {
    if (is.list(expr))
        expr = expr[[record$platform]]

    # get expression values from source name
    pd = Biobase::pData(expr)
    mapping = setNames(rownames(pd), pd$Source.Name)
    control = Biobase::exprs(expr)[,mapping[record$control]]
    perturbed = Biobase::exprs(expr)[,mapping[record$perturbed]]

    # build models
    mean_control= apply(control, 1, mean)
    sd_control = apply(control, 1, sd)
    fc = mean_control / apply(perturbed, 1, mean)
    model = loess(sd_control ~ mean_control)
    zscores = fc / predict(model, mean_control)

    # create index column
    index = record
    index$control = paste(index$control, collapse=";")
    index$perturbed = paste(index$perturbed, collapse=";")

    list(zscores=zscores, index=index)
}

records = lapply(yaml, function(r) record2zscore(r, expr))
index = lapply(records, function(x) x$index) %>% do.call(rbind, .)
zscores = lapply(records, function(x) x$zscores) %>% do.call(cbind, .)
result = list(index=index, zscores=zscores)

# separate index file w/ metadata derived from yaml [preferred?]
save(result, file=ZFILE)
