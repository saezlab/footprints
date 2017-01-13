library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
lincs = import('data/lincs')

GENESETS = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/reactome.RData"
INDEX = commandArgs(TRUE)[2] %or% "../../util/lincs/index.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "pathifier.RData"

row2scores = function(i, index, exps, sets) {
    library(dplyr)
    df = import('data_frame')
    lincs = import('data/lincs')
    pathifier = import_package('pathifier')

    row = index[i,]
    sign = row$sign
    ptb = df$subset(exps, row)$distil_id

    row$pathway = "control"
    row$pert_id = "DMSO"
    row$pert_dose = NULL
    row$sign = "0"
    ctl = df$subset(exps, row)$distil_id
    if (length(ctl) > 100)
        ctl = sample(ctl, 100, replace=FALSE)

    expr = lincs$get_z(cid=c(ctl,ptb), rid=lincs$projected, map_genes="hgnc_symbol")
    colnames(expr) = c(rep("ctl", length(ctl)), rep("ptb", length(ptb)))

    result = pathifier$quantify_pathways_deregulation(
        data = expr,
        allgenes = rownames(expr),
        syms = sets,
        pathwaynames = names(sets),
        normals = colnames(expr) == "ctl",
# default values
#        attempts = 100, # maybe set this higher and see if fewer NAs
        min_exp = -Inf,
#        min_std = 0.4
    )

    do.call(cbind, lapply(result$scores, c)) %>%
        ar$map(along=1, subsets=colnames(expr), mean) %>%
        ar$map(along=1, function(x) x['ptb']-x['ctl'])
}

# load model vectors and experiment index
sets = io$load(GENESETS)
exps = io$load(INDEX)

index = exps %>%
    select(pathway, cell_id, pert_id, pert_dose, pert_time, sign) %>%
    filter(pathway != "control") %>%
    distinct()

scores = clustermq::Q(row2scores, i=seq_len(nrow(index)),
                      const=list(index=index, exps=exps, sets=sets),
                      memory=10240, n_jobs=50) %>%
    setNames(seq_len(nrow(index))) %>%
    ar$stack(along=1) %>%
    ar$map(along=1, scale)

#errors = sapply(result, function(r) class(r) == "try-error")
#if (any(errors)) {
#    print(result[errors])
#    result[errors] = NA
#}

save(scores, index, file=OUTFILE)
