b = import('base')
io = import('io')
gdsc = import('data/gdsc')
gsea = import('../../util/gsea')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "go.RData"

# load gene lists for pathways
genelist = io$load(INFILE)

# get index, expr data for test set
speed = io$load('../../data/expr.RData')
keep = sapply(speed$records, function(x) identical(x$exclusion, "test-set"))
index = speed$records[keep]
expr = speed$expr[keep]

# perform GSEA
scores = mapply(gsea$runGSEA, expr=expr,
    MoreArgs=list(sigs=genelist, transform.normal=TRUE))

save(scores, index, file=OUTFILE)
