library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "scores_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "linear.pdf"

# load scores
scores = io$load(INFILE)

# load sanger data
Ys = gdsc$getDrugResponse('IC50s') # or AUC
tissues = gdsc$getTissues(minN=5)
ar$intersect(scores, tissues, Ys, along=1)

# save pdf w/ pan-cancer & tissue specific
pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

st$lm(Ys ~ tissues + scores) %>%
    filter(term == "scores") %>%
    select(-term, -tissues) %>%
    mutate(p.adj = p.adjust(p.value, method="fdr"),
           label = paste(Ys, scores, sep=":")) %>%
    plt$color$p_effect(pvalue="p.adj", effect="estimate", dir=-1) %>%
    plt$volcano(base.size=0.2) %>%
    print()

##TODO: filter by screening range
#Yf = sg$filterDrugResponse(Ys, tissues, top=0.1, abs=0, delta=2) # sensitive cell lines
#assocs.tissue = st$lm(Yf, scores, subsets=tissues, p.adjust="fdr", stack=T)
#print(plt$volcano(assocs.tissue, top=40, p=0.2, log='y', base.size=2))
