library(dplyr)
library(ArrayExpress)
b = import('base')
io = import('io')
ma = import('process/microarray')

dset = ArrayExpress('E-GEOD-15490') %>%
    ma$normalize() %>%
    ma$annotate('hgnc_symbol')

samples = io$read_table('samples.txt', header=FALSE)
colnames(samples) = c('id', 'name')
response = io$read_table('response.txt', header=TRUE) %>%
    mutate(patient = `Sample ID`)

expr = exprs(dset)
meta = pData(dset) %>%
    add_rownames('id') %>%
    transmute(id = sub(".CEL", "", id),
              treatment = Characteristics.treatment.) %>%
    inner_join(samples, by="id") %>%
    mutate(treated = grepl("after", treatment),
           treatment = b$grep("(FC|RFC|Rituximab)", treatment),
           rep = as.integer(sub(".*_rep", "", name)),
           patient = paste0("CLL", rep + ifelse(treatment != "FC", 10, 0))) %>%
    inner_join(response, by="patient") %>%
    dplyr::select(-rep, -name, -Treatment, -`Sample ID`)

save(expr, meta, file="data.RData")
