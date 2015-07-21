b = import('base')
io = import('io')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../genesets/reactome.RData"
EXPR = commandArgs(TRUE)[2] %or% "../../expr_cluster/corrected_expr.h5"
OUTFILE = commandArgs(TRUE)[3] %or% "pathifier.RData"

tissue2scores = function(tissue, EXPR) {
    io = import('io')
    pathifier = import_package('pathifier')

    tissues = io$h5load(EXPR, "/tissue")
    tumors = t(io$h5load(EXPR, "/expr", index=which(tissues == tissue)))
    normals = t(io$h5load(EXPR, "/expr", index=which(tissues == paste0(tissue, "_N"))))

    pathifier$quantify_pathways_deregulation(
        data = cbind(tumors, normals),
        allgenes = rownames(tumors),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = c(rep(FALSE, ncol(tumors)), rep(TRUE, ncol(normals))),
        attempts = 100,
        min_exp = -Inf,
        min_std = 0.4
    )
}

# load pathway gene sets
genesets = io$load(INFILE)

# get all tissues which have a normal
tissues = io$h5load(EXPR, "/tissue")
tissues = sub("_N", "", unique(tissues[grepl("_N", tissues)]))

# run pathifier in jobs
result = hpc$Q(tissue2scores, tissue=tissues, more.args=list(EXPR=EXPR), memory=4096)

# save results
save(result, file=OUTFILE)
