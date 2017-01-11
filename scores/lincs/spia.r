library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
lincs = import('data/lincs')

INDEX = commandArgs(TRUE)[1] %or% "../../util/lincs/index.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.RData"

row2scores = function(i) {
    row = index[i,]
    sign = row$sign
    ptb = df$subset(exps, row)$distil_id

    row$pathway = "control"
    row$pert_id = "DMSO"
    row$pert_dose = NULL
    row$sign = "0"
    ctl = df$subset(exps, row)$distil_id

### SPIA
    spia = import('../../util/spia')
    b = import('base/operators')

    print(rec)

    rownames(exp) = lookup$entrezgene[match(rownames(exp), lookup$hgnc_symbol)]
    exp = limma::avereps(exp[!is.na(rownames(exp)),])

    result = spia$spia(exp[,rec$perturbed], exp[,rec$control], pathids=spia$speed2kegg)

#    # if per_sample=TRUE
#    colnames(result) = spia$kegg2speed[colnames(result)]
#    result

    setNames(result, spia$kegg2speed[names(result)])
### SPIA END

    expr_ctl = expr[ctl,,drop=FALSE]
    expr_ptb = expr[ptb,,drop=FALSE]
    if (sign == "+")
        colMeans(expr_ptb) - colMeans(expr_ctl)
    else
        colMeans(expr_ctl) - colMeans(expr_ptb)
}

# load model vectors and experiment index
exps = io$load(INDEX)

#HGNC -> entrez gene lookup
lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
    biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"), mart=.)


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
#scores = ar$stack(result, along=1)

save(scores, index, file=OUTFILE)
