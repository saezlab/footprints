b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')
gsea = import('../../util/gsea')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathways_mapped/gsea_go.RData"

# only filter when we didn't select manually
if (grepl("mapped", OUTFILE)) {
    MIN_GENES = 1
    MAX_GENES = Inf
} else {
    MIN_GENES = 5
    MAX_GENES = 500
}

# load gene expression data, make sure same genes and drop duplicates
expr = lapply(tcga$tissues(), tcga$rna_seq) %>%
    ar$stack(along=2)
genelist = io$load(INFILE) %>%
    gsea$filter_genesets(rownames(expr), MIN_GENES, MAX_GENES)

# perform GSEA for each sample and signature
result = hpc$Q(gsea$runGSEA, sigs=genelist,
               const = list(expr=expr, transform.normal=TRUE),
               memory = 20480, job_size = 50)

# assemble results
result = setNames(result, names(genelist)) %>%
    ar$stack(along=2)

save(result, file=OUTFILE)
