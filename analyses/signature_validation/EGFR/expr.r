b = import('base')
ma = import('process/microarray')

expr = ArrayExpress::ArrayExpress("E-GEOD-6784") %>%
    ma$qc() %>%
    ma$normalize() %>%
    ma$annotate(summarize="hgnc_symbol")

if (length(expr) == 1)
    expr = expr[[1]]

save(expr, file="E-GEOD-6784.RData")
