b = import('base')
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')
tcga = import('data/tcga')
gsea = import('../../genesets/gsea')

INFILE = commandArgs(TRUE)[1] %or% "../../genesets/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "go.RData"

# load gene list and expression
genelist = io$load(INFILE)

tt = tcga$tissues()
expr = lapply(tt, tcga$rna_seq)
for (i in seq_along(expr)) {
    stopifnot(rownames(expr[[1]]) == rownames(expr[[i]]))
    expr[[i]] = expr[[i]][,!duplicated(colnames(expr[[i]]))]
}

# perform GSEA and save result
result = lapply(expr, function(e)
    gsea$runGSEA(e, genelist, transform.normal=TRUE))
for (r in result)
    stopifnot(colnames(result[[1]]) == colnames(r))
result = do.call(rbind, result)
result = result[!duplicated(rownames(result)),]
save(result, file=OUTFILE)
