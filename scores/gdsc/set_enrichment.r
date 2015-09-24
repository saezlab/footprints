b = import('base')
io = import('io')
gdsc = import('data/gdsc')
gsea = import('../../util/gsea')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "gsea_go.RData"

# load gene list and expression
genelist = io$load(INFILE)
expr = gdsc$basal_expression()

# perform GSEA
result = gsea$runGSEA(expr, genelist, transform.normal=TRUE)

save(result, file=OUTFILE)
