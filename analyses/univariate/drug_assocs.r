library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

# load scores
scores = io$load(INFILE)

# load sanger data
Ys = gdsc$drug_response('IC50s') # or AUC
tissues = gdsc$tissues(minN=5)
ar$intersect(scores, tissues, Ys, along=1)

# save pdf w/ pan-cancer & tissue specific
pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

# tissues as covariate
st$lm(Ys ~ tissues + scores) %>%
    filter(term == "scores") %>%
    select(-term, -tissues) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"),
           label = paste(Ys, scores, sep=":")) %>%
    plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
    plt$volcano(base.size=0.2) %>%
    print()

# separate associations for each tissue
Yf = gdsc$drug_response('IC50s', min_tissue_measured=2)
ar$intersect(Yf, scores, tissues)

st$lm(Yf ~ scores, subsets=tissues) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"),
           label = paste(subset, Yf, scores, sep=":")) %>%
    ungroup() %>%
    plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
    plt$volcano(p=0.2) %>%
    print()
