library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
lincs = import('data/lincs')

GENESETS = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/reactome.RData"
INDEX = commandArgs(TRUE)[2] %or% "../../util/lincs/index.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "speed_linear.RData"

row2scores = function(i) {
    row = index[i,]
    sign = row$sign
    ptb = df$subset(exps, row)$distil_id

    row$pathway = "control"
    row$pert_id = "DMSO"
    row$pert_dose = NULL
    row$sign = "0"
    ctl = df$subset(exps, row)$distil_id


### PATHIFIER
    data = exp[,c(rec$control, rec$perturbed)]
    colnames(data) = c(rep("control", length(rec$control)),
                       rep("perturbed", length(rec$perturbed)))

    result = pathifier$quantify_pathways_deregulation(
        data = data,
        allgenes = rownames(data),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = colnames(data) == "control",
# default values
#        attempts = 100, # maybe set this higher and see if fewer NAs
#        min_exp = 4,
#        min_std = 0.4
    )

    do.call(cbind, lapply(result$scores, c)) %>%
        ar$map(along=1, subsets=colnames(data), mean) %>%
        ar$map(along=1, function(x) x['perturbed']-x['control'])
### PATHIFIER END


    expr_ctl = expr[ctl,,drop=FALSE]
    expr_ptb = expr[ptb,,drop=FALSE]
    if (sign == "+")
        colMeans(expr_ptb) - colMeans(expr_ctl)
    else
        colMeans(expr_ctl) - colMeans(expr_ptb)
}

# load model vectors and experiment index
sets = io$load(GENESETS)
exps = io$load(INDEX)

index = exps %>%
    select(pathway, cell_id, pert_id, pert_dose, pert_time, sign) %>%
    filter(pathway != "control") %>%
    distinct()

scores = pbapply::pblapply(seq_len(nrow(index)), row2scores) %>%
    setNames(seq_len(nrow(index))) %>%
    ar$stack(along=1) %>%
    ar$map(along=1, scale)

#errors = sapply(result, function(r) class(r) == "try-error")
#if (any(errors)) {
#    print(result[errors])
#    result[errors] = NA
#}

save(scores, index, file=OUTFILE)
