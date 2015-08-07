b = import('base')
io = import('io')
gdsc = import('data/gdsc')
gsea = import('../../util/gsea')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "go.RData"

# load gene list and expression
genelist = io$load(INFILE)
speed = io$load('../../data/dscores.RData')

expr = speed$scores
index = speed$index %>%
    dplyr::select(-control, -perturbed)

# perform GSEA
scores = gsea$runGSEA(expr, genelist, transform.normal=TRUE)

save(scores, index, file=OUTFILE)
