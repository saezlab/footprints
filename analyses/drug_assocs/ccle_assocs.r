library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
ccle = import('data/ccle')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/ccle/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_matrix.RData"

# load required data
scores = io$load(INFILE)
drug = ccle$drug_response()
tissues = ccle$tissues(minN=5)
ar$intersect(scores, tissues, drug, along=1)

# tissues as covariate
pan = st$lm(drug ~ tissues + scores) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

# tissues as subsets
tissue = st$lm(drug ~ scores, subsets=tissues) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

save(pan, tissue, file=OUTFILE)
