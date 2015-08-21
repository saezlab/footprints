library(dplyr)
io = import('io')

mut = io$read_table("mutations_annotated_pathwayactivities.txt",
        na.strings=c("NA", "None", ""), header=TRUE) %>%
    select(HGNC = GENE_NAME,
           barcode = Tumor_Sample_Barcode,
           consequence = geneVeredict,
           H2O2, `IL-1`, `JAK-STAT`, MAPK, EGFR, PI3K, TGFb,
           TNFa, VEGF, Wnt, Insulin, Hypoxia, RAR, p53, Estrogen,
           Trail, notch, NFkB, PPAR, SHH) %>%
    filter(consequence %in% c("gainOfFunction", "lossOfFunction")) %>%
    tidyr::gather(path, mutated, -HGNC, -barcode, -consequence) %>%
    filter(!is.na(mutated)) %>%
    filter(!(HGNC=="TP53" && path!="p53")) %>%
    mutate(consequence = sub("OfFunction", "", consequence),
           barcode = substr(barcode, 1, 16)) %>%
    distinct()

print(sort(table(mut$HGNC)))
print(sort(table(mut$path)))

save(mut, file="mut_filtered.RData")
