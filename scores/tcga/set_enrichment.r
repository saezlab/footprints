b = import('base')
io = import('io')
gdsc = import('data/gdsc')
tcga = import('data/tcga')
gsea = import('../../genesets/gsea')

INFILE = commandArgs(TRUE)[1] %or% "../../genesets/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "go.RData"

# load gene list and expression
genelist = io$load(INFILE)

tt = tcga$tissues()
expr = lapply(tt, tcga$rna_seq)
for (e in expr)
    stopifnot(rownames(expr[[1]]) == rownames(e))
expr = do.call(cbind, expr)
expr = expr[,!duplicated(colnames(expr))]

# perform GSEA
result = gsea$runGSEA(expr, genelist, transform.normal=TRUE)

save(result, file=OUTFILE)
