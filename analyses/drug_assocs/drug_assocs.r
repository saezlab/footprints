library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "gdsc/speed_linear.pdf"

# load scores
scores = io$load(INFILE)

# load sanger data
if (grepl("gdsc", INFILE)) {
    gdsc = import('data/gdsc')
    Ys = gdsc$drug_response('IC50s') # or AUC
    Yf = gdsc$drug_response('IC50s', min_tissue_measured=2)
    tissues = gdsc$tissues(minN=5)
    ar$intersect(scores, tissues, Ys, Yf, along=1)
} else if (grepl("ccle", INFILE)) {
    ccle = import('data/ccle')
    Ys = ccle$drug_response()
    Yf = Ys
    Yf[Ys == 8] = NA
    tissues = ccle$tissues() #TODO:
    ar$intersect(scores, tissues, Ys, Yf, along=1)
}

# save pdf w/ pan-cancer & tissue specific
pdf(OUTFILE, paper="a4r", width=26, height=20)

# tissues as covariate
assocs.pan = st$lm(Ys ~ tissues + scores) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

# volcano plot for pan-cancer
assocs.pan %>%
    mutate(label = paste(Ys, scores, sep=":")) %>%
    plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
    plt$volcano(base.size=0.2) %>%
    print()

# matrix plot for pan-cancer
assocs.pan %>%
    mutate(lp = -log(adj.p),
           label = ifelse(adj.p < 1e-2, '*', ''),
           estimate = ifelse(adj.p < 0.1, estimate, NA)) %>%
    plt$cluster(lp ~ scores + Ys, size=c(Inf,20)) %>%
    plt$matrix(estimate ~ scores + Ys)

# tissues as subsets
assocs.tissue = st$lm(Yf ~ scores, subsets=tissues) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

# volcano plot for tissue subsets
assocs.tissue %>%
    mutate(label = paste(subset, Yf, scores, sep=":")) %>%
    plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
    plt$volcano(p=0.2) %>%
    print()

## matrix plots for tissue subsets # this doesn't really work, not many assocs
#assocs.tissue %>%
#    mutate(lp = -log(adj.p),
#           id = paste(scores, subset, sep=":"),
#           label = ifelse(adj.p < 0.05, '*', ''),
#           estimate = ifelse(adj.p < 0.2, estimate, NA)) %>%
#    plt$cluster(lp ~ id + Yf, size=c(30,40)) %>%
#    plt$matrix(estimate ~ id + Yf)

#TODO: linear fits of significant pan-can assocs

dev.off()
