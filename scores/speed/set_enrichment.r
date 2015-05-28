b = import('base')
io = import('io')
gdsc = import('data/gdsc')
icgc = import('data/icgc')
gsea = import('../../genesets/gsea')

INFILE = commandArgs(TRUE)[1] %or% "../../genesets/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "go.RData"

# load gene list and expression
genelist = io$load(INFILE)
speed = io$load('../../data/dscores.RData')
expr = speed$scores

# perform GSEA
result = gsea$runGSEA(expr, genelist, transform.normal=TRUE)

save(result, file=OUTFILE)
