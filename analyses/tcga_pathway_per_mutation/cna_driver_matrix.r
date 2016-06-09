library(dplyr)
io = import('io')
ar = import('array')
tcga = import('data/tcga')
gdsc = import('data/gdsc')

drivers = unique(gdsc$drivers()$HGNC)
cnat = tcga$cna_thresholded() %>%
    filter(hgnc %in% drivers) %>%
    mutate(barcode = substr(barcode, 1, 16),
           hgnc = ifelse(gistic > 0, paste0(hgnc, "_amp"), paste0(hgnc, "_del")),
           gistic = abs(gistic))

#FIXME: this loses 500/10k cases because no CNAs in drivers
cna_mat = ar$construct(gistic ~ barcode + hgnc, cnat, fill=0, fun.aggregate=mean)

save(cna_mat, file="cna_driver_matrix.RData")
