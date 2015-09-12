# todo: move somewhere else(?)
records2pathway = function(recs, expr, genesets) {
    library(dplyr)
    import('base/operators')
    ar = import('array')
    pathifier = import_package('pathifier')

    record2pathway = function(rec) {
        print(rec)

        data = expr[,c(rec$control, rec$perturbed)]
        colnames(data) = c(rep("control", length(rec$control)),
                           rep("perturbed", length(rec$perturbed)))

        result = pathifier$quantify_pathways_deregulation(
            data = data,
            allgenes = rownames(data),
            syms = genesets,
            pathwaynames = names(genesets),
            normals = colnames(data) == "control",
            attempts = 100,
            min_exp = -Inf,
            min_std = 0.4
        ) %catch% NA

        if (is.na(result))
            return(NA)

        do.call(cbind, lapply(result$scores, c)) %>%
            ar$map(along=1, subsets=colnames(data), mean) %>%
            ar$map(along=1, function(x) x['perturbed']-x['control'])
    }

    lapply(recs, record2pathway)
}

if (is.null(module_name())) {
    library(dplyr)
    io = import('io')
    hpc = import('hpc')

    OUTFILE = "pathifier_scores.RData"

    genesets = io$load("../util/genesets/reactome.RData")
    data = io$load('expr.RData')
    expr = data$expr
    records = data$records

    scores = hpc$Q(records2pathway,
                   recs = records, expr = expr,
                   const = list(genesets=genesets),
                   memory=4096, n_jobs=50) %>%
        unlist(recursive=FALSE)

    records = unlist(records, recursive=FALSE)
    stopifnot(length(records) == length(scores))

    save(scores, records, file=OUTFILE)
}
