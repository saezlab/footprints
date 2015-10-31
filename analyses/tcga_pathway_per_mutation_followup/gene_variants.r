b = import('base')
io = import('io')
tcga = import('data/tcga')
cbp = import('data/cbioportal')

INFILE = commandArgs(TRUE)[1] %or% "../tcga_pathway_per_mutation/mutations_annotated_pathwayactivities_v3_mikeformat.txt"
OUTFILE = commandArgs(TRUE)[2] %or% "gene_variants.RData"

genes = io$read_table(INFILE, header=TRUE) %>%
    group_by("GENE_NAME") %>%
    filter(n() > 10) %>%
    ungroup() %>%
    select(GENE_NAME) %>%
    unlist() %>%
    unique()

# for all tcga studies
studies = cbp$studies(tcga=TRUE)

# get SNVs
variants = cbp$variants(studies, genes) %>%
    transmute(id = id,
              study = tcga$barcode2study(id),
              hgnc = hgnc,
              variant = amino_acid_change) %>%
    distinct() %>%
    na.omit()

# and CNAs
lookup = c('-2' = "homo_del",
           '-1' = "hemi_del",
           '0' = NA,
           '1' = "low_amp",
           '2' = "high_amp")

cna = cbp$cna(studies, genes, only_altered=TRUE) %>%
    tidyr::gather(cna, value, -id) %>%
    mutate(cna = as.character(cna),
           value = unlist(lookup[as.character(value)])) %>%
    transmute(id = id,
              study = tcga$barcode2study(id),
              hgnc = cna,
              variant = value) %>%
    distinct() %>%
    na.omit()

# save to object
save(variants, cna, file=OUTFILE)
