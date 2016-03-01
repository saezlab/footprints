b = import('base')
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')
gsea = import('../../util/gsea')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathways_mapped/gsea_go.RData"
MIN_GENES = 5
MAX_GENES = 500

# load gene list and expression
genelist = io$load(INFILE)
expr = gdsc$basal_expression()

# filter gene list by number of genes
num_overlap = sapply(genelist, function(x) length(intersect(rownames(expr), x)))
discard = num_overlap < MIN_GENES | num_overlap > MAX_GENES
if (any(discard)) {
    warning("Discarding the following sets: ", paste(names(genelist)[discard], collapse=", "))
    genelist = genelist[!discard]
}

# perform GSEA
result = hpc$Q(gsea$runGSEA, sigs=genelist,
               const = list(expr=expr, transform.normal=TRUE),
               memory = 4096, job_size = 50)

result = setNames(result, names(genelist)) %>%
    ar$stack(along=2)

save(result, file=OUTFILE)
