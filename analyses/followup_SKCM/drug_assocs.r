library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_matrix.pdf"

# response data for 52 cell lines
# >=10 measured data points: 76 drugs
# >=20 measured data points: 45 drugs
scores = io$load(INFILE)
tissues = gdsc$tissues("SKCM")
Yf = gdsc$drug_response('IC50s', min_tissue_measured=10)
ar$intersect(scores, tissues, Yf, along=1)
Yf = Yf[,ar$map(Yf, along=1, function(x) !all(is.na(x)))]

# tissues as subsets
assocs.tissue = st$lm(Yf ~ scores) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

# volcano plot for tissue subsets
pdf(OUTFILE, paper="a4r", width=26, height=20)
assocs.tissue %>%
    mutate(label = paste(Yf, scores, sep=":")) %>%
    plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
    plt$volcano(p=0.2) %>%
    print()
dev.off()
