b = import('base')
io = import('io')
gdsc = import('data/gdsc')
gsea = import('../../util/gsea')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathways_mapped/gsea_go.RData"

# load gene list and expression
genelist = io$load(INFILE)
expr = gdsc$basal_expression()

# perform GSEA
result = hpc$Q(gsea$runGSEA, sigs=genelist,
               const = list(expr=expr, transform.normal=TRUE),
               memory = 2048, job_size = 50)

result = setNames(result, names(genelist)) %>%
    ar$stack(along=2)

save(result, file=OUTFILE)
