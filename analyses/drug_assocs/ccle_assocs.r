library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
ccle = import('data/ccle')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.RData"

# load required data
scores = io$load(INFILE)
Ys = ccle$drug_response()
Yf = Ys
Yf[Ys == 8] = NA
tissues = ccle$tissues() #TODO:
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
