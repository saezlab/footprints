import('base/operators')
io = import('io')
ar = import('array')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/reactome.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "pathifier.RData"

tissue2scores = function(tissue, genesets) {
    io = import('io')
    tcga = import('data/tcga')
    pathifier = import_package('pathifier')

    expr = tcga$rna_seq(tissue)
    is_normal = grepl("[Nn]ormal", tcga$barcode2index(colnames(expr))$Sample.Definition)

    result = pathifier$quantify_pathways_deregulation(
        data = expr,
        allgenes = rownames(expr),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = is_normal,
        attempts = 100,
        min_exp = -Inf,
        min_std = 0.4
    )

    re = do.call(cbind, lapply(result$scores, c))
    rownames(re) = colnames(data)
    re
}

# load pathway gene sets
genesets = io$load(INFILE)
tissues = c("BLCA", "BRCA", "CESC", "COREAD", "ESCA", "HNSC",
            "KIRC", "LIHC", "LUAD", "LUSC", "PAAD")

# run pathifier in jobs
result = hpc$Q(tissue2scores, tissue=tissues,
    const=list(genesets=genesets), memory=8192, n_jobs=length(tissues))

result = ar$stack(result, along=1)

# save results
save(result, file=OUTFILE)
