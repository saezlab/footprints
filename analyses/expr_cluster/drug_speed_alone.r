library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/merge/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "drug_speed_alone.pdf"

# load scores
scores = io$load(INFILE)
scores = scores[grepl("^[0-9]+$", rownames(scores)),]

# load sanger data
# separate associations for each tissue
Yf = gdsc$drug_response('IC50s', min_tissue_measured=2)
tissues = gdsc$tissues()

ar$intersect(scores, Yf, tissues)

# save pdf w/ pan-cancer & tissue specific
pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

for (tissue in unique(tissues)) {
    message(tissue)
    score = scores[tissues == tissue,]
    resp = Yf
    ar$intersect(score, resp, along=1)

    # tissues as subsets
    assocs.tissue = st$lm(resp ~ score) %>%
        filter(term == "score") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               label = paste(resp, score, sep=":")) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
        plt$volcano(p=0.2) + ggtitle(tissue)
    print(assocs.tissue)
}
