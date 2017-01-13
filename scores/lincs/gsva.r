library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
lincs = import('data/lincs')

GENESETS = commandArgs(TRUE)[1] %or% "../../util/genesets/mapped/go.RData"
INDEX = commandArgs(TRUE)[2] %or% "../../util/lincs/index.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "gsva_go.RData"

row2scores = function(i) {
    row = index[i,]
    sign = row$sign
    ptb = df$subset(exps, row)$distil_id

    row$pathway = "control"
    row$pert_id = "DMSO"
    row$pert_dose = NULL
    row$sign = "0"
    ctl = df$subset(exps, row)$distil_id

    gsva_ctl = gsva[ctl,,drop=FALSE]
    gsva_ptb = gsva[ptb,,drop=FALSE]
    colMeans(gsva_ptb) - colMeans(gsva_ctl)
}

# load model vectors and experiment index
sets = io$load(GENESETS)
exps = io$load(INDEX)
expr = lincs$get_z(exps$distil_id, rid=lincs$projected, map_genes="hgnc_symbol")

# calculate scores for all samples
gsva_fun = function(set, sets, expr)
    GSVA::gsva(expr=expr, gset.idx.list=sets[set], parallel.sz=1)$es.obs
gsva = clustermq::Q(gsva_fun, set=names(sets), const=list(sets=sets, expr=expr),
                    memory=10240, job_size=1) %>% ar$stack(along=2) %>% t()

index = exps %>%
    select(pathway, cell_id, pert_id, pert_dose, pert_time, sign) %>%
    filter(pathway != "control") %>%
    distinct()

scores = pbapply::pblapply(seq_len(nrow(index)), row2scores) %>%
    setNames(seq_len(nrow(index))) %>%
    ar$stack(along=1) %>%
    ar$map(along=1, scale)

save(scores, index, file=OUTFILE)
