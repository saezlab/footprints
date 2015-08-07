tissue2scores = function(tissue, genesets, EXPR) {
    io = import('io')
    pathifier = import_package('pathifier')

    tissues = io$h5load(EXPR, "/tissue")
    tumors = t(io$h5load(EXPR, "/expr", index=which(tissues == tissue)))
    normals = t(io$h5load(EXPR, "/expr", index=which(tissues == paste0(tissue, "_N"))))
    data = cbind(tumors, normals)

    result = pathifier$quantify_pathways_deregulation(
        data = data,
        allgenes = rownames(tumors),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = c(rep(FALSE, ncol(tumors)), rep(TRUE, ncol(normals))),
        attempts = 100,
        min_exp = -Inf,
        min_std = 0.4
    )

    re = do.call(cbind, lapply(result$scores, c))
    rownames(re) = colnames(data)
    re
}

if (is.null(module_name())) {
    b = import('base')
    io = import('io')
    ar = import('array')
    hpc = import('hpc')

    INFILE = commandArgs(TRUE)[1] %or% "../../genesets/reactome.RData"
    EXPR = commandArgs(TRUE)[2] %or% "../../expr_cluster/corrected_expr.h5"
    OUTFILE = commandArgs(TRUE)[3] %or% "pathifier.RData"

    # load pathway gene sets
    genesets = io$load(INFILE)

    # get all tissues which have a normal
    tissues = io$h5load(EXPR, "/tissue")
    tissues = sub("_N", "", unique(tissues[grepl("_N", tissues)]))

    # run pathifier in jobs
    result = hpc$Q(tissue2scores, tissue=tissues,
        more.args=list(EXPR=EXPR, genesets=genesets), memory=8192)

    result = ar$stack(result, along=1)

    # save results
    save(result, file=OUTFILE)
}
