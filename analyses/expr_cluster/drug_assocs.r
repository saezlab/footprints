library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "../../expr_cluster/speed_cluster.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "drug_speed.pdf"

# load scores
clusters = io$load(INFILE)
clusters = clusters[lapply(clusters, length) != 0]
scores = lapply(clusters, function(cl) {
    re = cl %>%
        lapply(function(x) x[grepl("^[0-9]+$", rownames(x)),]) %>%
        ar$stack(along=2)
    re[,colSums(re)>0]
})

# load sanger data
# separate associations for each tissue
Yf = gdsc$drug_response('IC50s', min_tissue_measured=2)
tissues = gdsc$tissues()

# save pdf w/ pan-cancer & tissue specific
pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

for (tissue in names(scores)) {
    message(tissue)
    score = scores[[tissue]]
    resp = Yf
    ar$intersect(score, resp, along=1)

    # require at least 2 cell lines in a cluster
    score = score[,colSums(score) > 1]

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
