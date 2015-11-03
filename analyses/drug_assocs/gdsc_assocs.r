library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.RData"

# load required data
scores = io$load(INFILE)
Ys = gdsc$drug_response('IC50s') # or AUC
Yf = gdsc$drug_response('IC50s', min_tissue_measured=5)
tissues = gdsc$tissues(minN=15)
ar$intersect(scores, tissues, Ys, Yf, along=1)

# tissues as covariate
assocs.pan = st$lm(Ys ~ tissues + scores) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

# tissues as subsets
assocs.tissue = st$lm(Yf ~ scores, subsets=tissues) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

save(assocs.pan, assocs.tissue, file=OUTFILE)
