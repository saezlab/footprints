library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
df = import('data_frame')
plt = import('plot')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "SCLC_assocs.pdf"

scores = io$load(INFILE)
tissues = gdsc$tissues("SCLC")
Yf = gdsc$drug_response('IC50s', min_tissue_measured=5, drop=TRUE)
mut = gdsc$mutated_genes(intogen=TRUE, tissue="SCLC", drop=TRUE) + 0
ar$intersect(scores, tissues, Yf, along=1)

# tissues as subsets
assocs.tissue = st$lm(Yf ~ scores) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

assocs.tissue %>%
    mutate(label = paste(Yf, scores, sep=":")) %>%
    plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
    plt$volcano(p=0.2)

dd = df$assemble(drug=Yf[,'PD-0325901'] , MAPK=scores[,'MAPK'] , TNFa=scores[,'TNFa'])
st$lm(drug ~ MAPK + TNFa, data=dd) # MAPK: p=0.046, TNFa:7e-4
