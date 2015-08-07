b = import('base')
io = import('io')
lincs = import('data/lincs')
gsea = import('../../util/gsea')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/go.RData"
INDEX = commandArgs(TRUE)[2] %or% "../../util/lincs_perturbation_qc/index.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "go.RData"

# load gene list and expression
genelist = io$load(INFILE)
index = unique(io$load(INDEX)$distil_id)
expr = lincs$get_z(cid=index, rid=lincs$projected, map.genes="hgnc_symbol")

# perform GSEA
scores = gsea$runGSEA(expr, genelist, transform.normal=TRUE)

#save(scores, index, file=OUTFILE)
save(scores, file=OUTFILE)
