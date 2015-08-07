b = import('base')
io = import('io')
ar = import('array')
spia = import('../../util/spia')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/reactome.RData"
EXPR = commandArgs(TRUE)[2] %or% "../../util/expr_cluster/corrected_expr.h5"
OUTFILE = commandArgs(TRUE)[3] %or% "spia.RData"

tissue2scores = function(tissue, EXPR, spia) {
    io = import('io')

    tissues = io$h5load(EXPR, "/tissue")
    tumors = t(io$h5load(EXPR, "/expr", index=which(tissues == tissue)))
    normals = t(io$h5load(EXPR, "/expr", index=which(tissues == paste0(tissue, "_N"))))

    re = spia$spia(tumors, normals, per_sample=TRUE, pathids=spia$speed2kegg, verbose=TRUE)
    setNames(re$score, re$name)
}

# load pathway gene sets and tissues
EXPR = tools:::file_path_as_absolute(EXPR)
genesets = io$load(INFILE)
tissues = io$h5load(EXPR, "/tissue")
tissues = sub("_N", "", unique(tissues[grepl("_N", tissues)]))

# run spia in jobs and save
result = hpc$Q(tissue2scores, tissue=tissues, more.args=list(EXPR=EXPR, spia=spia), memory=4096)
result = ar$stack(result, along=1)
save(result, file=OUTFILE)
