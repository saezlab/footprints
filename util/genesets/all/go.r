library(biomaRt)
library(modules)
library(dplyr)
b = import('base')

OUTFILE = commandArgs(TRUE)[1] %or% "go.RData"
MIN_GENES = 5
MAX_GENES = 100

# 319k x2 for all genes, 20k x2 for landmarks
mart = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

mapGO = biomaRt::getBM(attributes=c("hgnc_symbol", "go_id", "name_1006"), mart=mart)

stacked = mapGO %>%
    filter(go_id != "") %>%
    group_by(go_id) %>%
    filter(n()>=MIN_GENES & n()<=MAX_GENES) %>%
    mutate(go=paste(go_id,name_1006)) %>%
    select(hgnc_symbol,go_id)

genesets = unstack(stacked)

save(genesets, file=OUTFILE)
