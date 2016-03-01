b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')
hpc = import('hpc')

#' Calculate GSEA score for one tissue and gene set
#'
#' @param genelist   A gene set of list of gene sets
#' @param tissue     The TCGA tissue to subset expression to
#' @param expr_list  A list of TCGA gene expression
tissue2scores = function(genelist, tissue, expr_list) {
    gsea = import('../../util/gsea')
    e = expr[[tissue]]
    gsea$runGSEA(e, genelist, transform.normal=TRUE)
}

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathways_mapped/go.RData"
MIN_GENES = 5
MAX_GENES = 500

# load gene expression data, make sure same genes and drop duplicates
genelist = io$load(INFILE)
expr = lapply(tcga$tissues(), tcga$rna_seq)
for (i in seq_along(expr)) {
    stopifnot(rownames(expr[[1]]) == rownames(expr[[i]]))
    expr[[i]] = expr[[i]][,!duplicated(colnames(expr[[i]]))]
}

# filter gene list by number of genes
num_overlap = sapply(genelist, function(x) length(intersect(rownames(expr[[1]]), x)))
discard = num_overlap < MIN_GENES | num_overlap > MAX_GENES
if (any(discard)) {
    warning("Discarding the following sets: ", paste(names(genelist)[discard], collapse=", "))
    genelist = genelist[!discard]
}

# calc GSEA for each tissue and gene set pair
result = hpc$Q(tissue2scores, genelist=genelist, tissue=names(expr),
               const = list(expr_list=expr),
               memory = 10240, job_size = 50, expand_grid=TRUE)

# assemble results
for (r in result)
    stopifnot(colnames(result[[1]]) == colnames(r))
result = do.call(rbind, result)
result = result[!duplicated(rownames(result)),]
save(result, file=OUTFILE)
