b = import('base')
io = import('io')
gdsc = import('data/gdsc')
icgc = import('data/icgc')
gsea = import('../../genesets/gsea')

INFILE = commandArgs(TRUE)[1] %or% "../../genesets/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "go.RData"

# load gene list and expression
genelist = io$load(INFILE)
expr = icgc$rna_seq(voom=TRUE, map_ids="icgc_specimen_id")
expr = expr[,!duplicated(colnames(expr))]

# perform GSEA
result = gsea$runGSEA(expr, genelist, transform.normal=TRUE)

save(result, file=OUTFILE)
