b = import('base')
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')
genesets = import('../../util/genesets')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathways_mapped/gsva_go.RData"

# only filter when we didn't select manually
if (grepl("mapped", OUTFILE)) {
    MIN_GENES = 1
    MAX_GENES = Inf
    job_size = 1
} else {
    MIN_GENES = 5
    MAX_GENES = 500
    job_size = 20
}

# load gene list and expression
expr = gdsc$basal_expression()
genesets = io$load(INFILE) %>%
    genesets$filter(rownames(expr), MIN_GENES, MAX_GENES)

#' Function to calculate GSVA score for a single signature
#'
#' @param set     A character describing the gene set (element in sigs)
#' @param expr    The gene expression matrix [genes x samples]
#' @param sigs    The list of signatures
#' @return        Result for GSVA(expr[,sample], sigs[set])
gsva = function(set, expr, sigs, ...) {
    GSVA::gsva(expr=expr, gset.idx.list=sigs[set], parallel.sz=0, ...)$es.obs
}

# perform GSEA for each sample and signature
result = hpc$Q(gsva, set = names(genesets),
               const = list(expr=expr, sigs=genesets),
               memory = 10240, job_size = job_size) %>%
    ar$stack(along=1) %>% t()

save(result, file=OUTFILE)
