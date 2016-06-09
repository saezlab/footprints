library(dplyr)
io = import('io')
ar = import('array')
tcga = import('data/tcga')
gdsc = import('data/gdsc')

drivers = unique(gdsc$drivers()$HGNC)
cnat = tcga$cna_thresholded() %>%
    mutate(barcode = substr(barcode, 1, 16)) %>%
    filter(hgnc %in% drivers)

#FIXME: this loses 500/10k cases because no CNAs in drivers
cna_mat = ar$construct(gistic ~ barcode + hgnc, cnat, fill=0, fun.aggregate=mean)

save(cna_mat, file="cna_driver_matrix.RData")
