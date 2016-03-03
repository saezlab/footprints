import('base/operators')
io = import('io')
ar = import('array')
tcga = import('data/tcga')
gsea = import('../../util/gsea')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/reactome.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "pathways_mapped/pathifier.RData"
MIN_GENES = 5
MAX_GENES = 500

tissue2scores = function(tissue, genesets, expr) {
    library(dplyr)
    io = import('io')
    tcga = import('data/tcga')
    pathifier = import_package('pathifier')

    index = tcga$barcode2index(colnames(expr)) %>%
        filter(Study.Abbreviation == tissue)
    my_expr = expr[,index$Bio.ID]
    is_normal = grepl("[Nn]ormal", index$Sample.Definition)

    result = pathifier$quantify_pathways_deregulation(
        data = my_expr,
        allgenes = rownames(my_expr),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = is_normal,
#        attempts = 100, # default settings
        min_exp = -Inf, #TODO: std=4, applicable to voom?
#        min_std = 0.4
    )

    re = do.call(cbind, lapply(result$scores, c))
    rownames(re) = colnames(my_expr)
    re
}

# load pathway gene sets
tissues = import('../../config')$tcga$tissues_with_normals
expr = lapply(tcga$tissues(), tcga$rna_seq) %>%
    ar$stack(along=2)
genesets = io$load(INFILE) %>%
    gsea$filter_genesets(rownames(expr), MIN_GENES, MAX_GENES)

# handle COAD and READ separately
if ("COADREAD" %in% tissues) {
    tissues = setdiff(tissues, "COADREAD")
    tissues = c(tissues, "COAD", "READ")
}

# make compatible to call with one set in above function
for (i in seq_along(genesets))
    genesets[[i]] = setNames(list(genesets[[i]]), names(genesets)[i])

# run pathifier in jobs
result = hpc$Q(tissue2scores, tissue=tissues, genesets=genesets,
               const = list(expr=expr), expand_grid=TRUE,
               memory=8192, job_size=5, fail_on_error=FALSE)

result[sapply(result, class) == "try-error"] = NA
result = ar$stack(result, along=2)

# save results
save(result, file=OUTFILE)
