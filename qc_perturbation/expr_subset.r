library(dplyr)
io = import('io')
#TODO: this should be saved in data/expr already in this format

OUTFILE = "expr_subset.RData"

data = io$load('../data/expr.RData')
expr = data$expr
records = unlist(data$records, recursive=FALSE, use.names=FALSE)
rec_names = sapply(records, function(x) x$id)

# use only records for the test set
keep = sapply(records, function(x) identical(x$exclusion, "test-set"))
records = setNames(records, rec_names)[keep]

record2expr = function(rec) {
    e = expr[[paste(rec$accession, rec$platform, sep=".")]]
    e[,c(rec$control, rec$perturbed)]
}
expr = setNames(lapply(records, record2expr), names(records))

save(records, expr, file=OUTFILE)
