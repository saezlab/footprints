library(dplyr)
io = import('io')
ma = import('process/microarray')

process_record = function(rec) {
    re = rec$expr
    message(re$accession, " ", re$platform)
    expr = ArrayExpress::ArrayExpress(re$accession, drop=FALSE)[[re$platform]] %>%
#        ma$qc() %>%
        ma$normalize() %>%
        ma$annotate(summarize="hgnc_symbol")
}

index = io$read_yaml("validation.yaml")
expr = lapply(index, process_record)
names(expr) = lapply(index, function(x) x$expr$accession)

save(expr, file="expr.RData")
