b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')
gdsc = import('data/gdsc')

studies = import('../../config')$tcga$tissues_with_normals
drivers = unique(gdsc$drivers()$HGNC)

# AAChange is not avail in eg. BRCA (and others)
mut = tcga$mutations(id_type="specimen", subset="primary") %>%
    filter(Study %in% studies &
           Variant_Classification != "Silent") %>%
    mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode)) %>%
    filter(Hugo_Symbol %in% drivers)

mut$has_mut = 1
mut_mat = ar$construct(has_mut ~ Tumor_Sample_Barcode + Hugo_Symbol,
                       data=mut, fun.aggregate = length) > 0

mut_mat = mut_mat + 0

save(mut_mat, file="mut_driver_matrix.RData")
