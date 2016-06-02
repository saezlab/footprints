b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')
gsea = import('../../util/gsea')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathways_mapped/gsva_go.RData"

# only filter when we didn't select manually
if (grepl("mapped", OUTFILE)) {
    MIN_GENES = 1
    MAX_GENES = Inf
} else {
    MIN_GENES = 5
    MAX_GENES = 500
}

# load gene list and expression
tissues = import('../../config')$tcga$tissues
expr = b$lnapply(tissues, tcga$rna_seq)
genesets = io$load(INFILE)

#' Function to calculate GSVA score for a single signature
#'
#' @param set     A character describing the gene set (element in sigs)
#' @param set     A character describing the tissue (element in expr)
#' @param expr    The gene expression matrix [genes x samples]
#' @param sigs    The list of signatures
#' @return        Result for GSVA(expr[,sample], sigs[set])
gsva = function(expr, sigs, ...) {
    GSVA::gsva(expr, gset.idx.list=sigs, parallel.sz=1, ...)$es.obs
}

# perform GSEA for each sample and signature
result = hpc$Q(gsva, expr=expr,
               const = list(sigs=genesets, min.sz=MIN_GENES, max.sz=MAX_GENES),
               memory = 10240, job_size = 1) %>%
    ar$stack(along=2) %>% t()

save(result, file=OUTFILE)
