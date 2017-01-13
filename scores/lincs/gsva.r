library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
lincs = import('data/lincs')

GENESETS = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/go.RData"
INDEX = commandArgs(TRUE)[2] %or% "../../util/lincs/index.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "gsva_go.RData"

row2scores = function(i, index, exps, sets) {
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
    gsva = GSVA::gsva(expr=expr, gset.idx.list=sets, parallel.sz=1)$es.obs

    gsva_ctl = gsva[ctl,,drop=FALSE]
    gsva_ptb = gsva[ptb,,drop=FALSE]
    colMeans(gsva_ptb) - colMeans(gsva_ctl)
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

save(scores, index, file=OUTFILE)
