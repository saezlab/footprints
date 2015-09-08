library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

# better: mapping yaml->zscore file
YAML = commandArgs(TRUE)[1] %or% list.files("new_index", "[0-9]+\\.yaml$",
                                 recursive=TRUE, full.names=TRUE)
EXPR = commandArgs(TRUE)[2] %or% list.files("normalized", "\\.RData", full.names=TRUE)
OUTFILE = commandArgs(TRUE)[3] %or% "./zscores_commonNULL.RData"

record2expr = function(record) {
    e = expr[[record$accession]]
    if (is.list(e))
        e = e[[record$platform]]

    pd = Biobase::pData(e)
    emat = Biobase::exprs(e)
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

    # create index column
    index = record
    index$control = paste(index$control, collapse=";")
    index$perturbed = paste(index$perturbed, collapse=";")

    list(index=index, control=control, perturbed=perturbed)
}

yaml = lapply(YAML, function(y) io$read_yaml(y, drop=FALSE)) %>%
    do.call(c, .)
expr = lapply(EXPR, function(e) io$load(e)) %>%
    setNames(b$grep("/([^/]+)\\.RData", EXPR))

all_combined = lapply(yaml, function(y) record2expr(y) %catch% NULL)
all_combined = all_combined[!sapply(all_combined, is.null)]

all_control = lapply(all_combined, function(x) x$control) %>%
    ar$stack(along=2) # this will fail w/ aggregation error, but doesn't matter in this case
mean_control = apply(all_control, 1, function(x) mean(x, na.rm=TRUE))
sd_control = apply(all_control, 1, function(x) sd(x, na.rm=TRUE))
model = loess(sd_control ~ mean_control)

record2zscore = function(record) {
    mean_perturbed = apply(record$perturbed, 1, function(x) mean(x, na.rm=TRUE))
    logFC = mean_perturbed - mean_control[names(mean_perturbed)]
    zscores = logFC / predict(model, mean_perturbed)
    list(index=record$index, dscores=logFC, zscores=zscores)   
}
records = lapply(all_combined, record2zscore)
names(records) = 1:length(records)
index = lapply(records, function(x) x$index) %>% bind_rows()
scores = lapply(records, function(x) x$zscores) %>% ar$stack(along=2)

save(index, scores, file=OUTFILE)
