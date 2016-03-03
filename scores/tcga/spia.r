b = import('base')
io = import('io')
ar = import('array')
hpc = import('hpc')
spia = import('../../util/spia')
tcga = import('data/tcga')

tissue2scores = function(tissue, pathids, expr) {
    io = import('io')
    tcga = import('data/tcga')
    spia = import('../../util/spia')

    index = tcga$barcode2index(colnames(expr)) %>%
        filter(Study.Abbreviation == tissue)
    my_expr = expr[,index$Bio.ID]
    is_normal = grepl("[Nn]ormal", index$Sample.Definition)

    tumors = my_expr[,!is_normal, drop=FALSE]
    normals = my_expr[,is_normal, drop=FALSE]

    spia$spia_per_sample(tumors, normals, pathids=spia$speed2kegg)
}

if (is.null(module_name())) {
    OUTFILE = commandArgs(TRUE)[1] %or% "pathways_mapped/spia.RData"
    FILTER = as.logical(commandArgs(TRUE)[2]) %or% TRUE

    # load pathway gene sets
    tissues = import('../../config')$tcga$tissues_with_normals
    expr = lapply(tcga$tissues(), tcga$rna_seq) %>%
        ar$stack(along=2) %>%
        spia$map_entrez()

    # handle COAD and READ separately
    if ("COADREAD" %in% tissues) {
        tissues = setdiff(tissues, "COADREAD")
        tissues = c(tissues, "COAD", "READ")
    }

    if (FILTER)
        genesets = spia$speed2kegg
    else
        genesets = spia$pathids("hsa")

    # make compatible to call with one set in above function
    for (i in seq_along(genesets))
        genesets[[i]] = setNames(list(genesets[[i]]), names(genesets)[i])

    # run spia in jobs and save
    result = hpc$Q(tissue2scores, tissue=tissues, pathids=genesets,
				   const = list(expr=expr), expand_grid=TRUE,
                   memory=8192, job_size=5, fail_on_error=FALSE) %>%
        setNames(tissues) %>%
        ar$stack(along=1)

    result[sapply(result, class) == "try-error"] = NA
    result = ar$stack(result, along=2)

    if (FILTER)
        colnames(result) = spia$kegg2speed[colnames(result)]

    save(result, file=OUTFILE)
}
