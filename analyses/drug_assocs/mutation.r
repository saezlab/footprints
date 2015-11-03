library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')

OUTFILE = commandArgs(TRUE)[1] %or% "mutation.pdf"

# load sanger data
scores = gdsc$mutated_genes(intogen=TRUE, drop=TRUE) + 0
Ys = gdsc$drug_response('IC50s') # or AUC
tissues = gdsc$tissues(minN=5)
ar$intersect(scores, tissues, Ys, along=1)

# tissues as covariate
assocs.pan = st$lm(Ys ~ tissues + scores) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

# separate associations for each tissue
tissue_assocs = function(tissue) {
    Yf = gdsc$drug_response('IC50s', min_tissue_measured=5)
    scores = gdsc$mutated_genes(intogen=TRUE, drop=TRUE) + 0
    tissues = gdsc$tissues(tissue)
    ar$intersect(Yf, scores, tissues)

    assocs.tissue = st$lm(Yf ~ scores) %>%
        filter(term == "scores") %>%
        select(-term) %>%
        mutate(subset = tissue,
               adj.p = p.adjust(p.value, method="fdr"))
}

assocs.tissue = unique(tissues) %>%
    lapply(tissue_assocs) %>%
    bind_rows()

save(assocs.pan, assocs.tissue, file=OUTFILE)
