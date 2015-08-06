tissue2scores = function(tissue, genesets, INDEX, EXPR) {
    io = import('io')
    pathifier = import_package('pathifier')

    # prepare data
    index = io$load(INDEX)
    expr = io$load(EXPR)
# tissues = pathway perturbation subsets
    tissues = index$pathway
# tumors = perturbed
    tumors = expr[,index$sign == "+" & index$pathway == tissue]
# normals = control
    normals = expr[,index$sign == "0" & index$pathway == tissue]
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
    INDEX = "../../data/lincs_perturbation_qc/index.RData"
    EXPR = commandArgs(TRUE)[2] %or% "../../data/lincs_perturbation_qc/expr.RData"
    OUTFILE = commandArgs(TRUE)[3] %or% "pathifier.RData"

    # load pathway gene sets
    genesets = io$load(INFILE)

#    # get all tissues which have a normal
#    tissues = io$h5load(EXPR, "/tissue")
#    tissues = sub("_N", "", unique(tissues[grepl("_N", tissues)]))
    tissues = unique(io$load(INDEX)$pathway)

    # run pathifier in jobs
    result = hpc$Q(tissue2scores, tissue=tissues,
        more.args=list(INDEX=INDEX, EXPR=EXPR, genesets=genesets), memory=4096)

    result = ar$stack(result, along=1)

    # save results
    save(result, file=OUTFILE)
}
