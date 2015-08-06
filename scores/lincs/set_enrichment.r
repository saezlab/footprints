b = import('base')
io = import('io')
gsea = import('../../genesets/gsea')

INFILE = commandArgs(TRUE)[1] %or% "../../genesets/go.RData"
EXPR = commandArgs(TRUE)[2] %or% "../../data/lincs_perturbation_qc/expr.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "go.RData"

# load gene list and expression
genelist = io$load(INFILE)
expr = io$load(EXPR)
#index = io$load('../../data/lincs_perturbation_qc/index.RData')

# perform GSEA
scores = gsea$runGSEA(expr, genelist, transform.normal=TRUE)

#save(scores, index, file=OUTFILE)
save(scores, file=OUTFILE)
