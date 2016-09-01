library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
st = import('stats')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "correlation.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "assocs.RData"

#' Calculate associations for one study/tissue and cor changes between mut/wt
#'
#' @param subs  The subset or study ("pan" or TCGA cancer type identifier)
#' @return      A data.frame with the mutation, delta correlation, and p.value
do_assoc = function(cor_row) {
    cur_mut = bem
    tissues = tcga$barcode2study(rownames(cur_mut))
    if (cor_row$tissue != "pan") {
        tissues = tissues[tissues == cor_row$tissue]
        ar$intersect(cur_mut, tissues, along=1)
    }

    mut2assoc = function(mut_hgnc) {
        is_wt = rownames(mut_hgnc)[!is.na(mut_hgnc[,1]) & mut_hgnc[,1] == 0]
        is_mut = rownames(mut_hgnc)[!is.na(mut_hgnc[,1]) & mut_hgnc[,1] == 1]

        pathways = c(cor_row$path1, cor_row$path2)
        cor_wt = scores[is_wt, pathways]
        cor_mut = scores[is_mut, pathways]

        # do correlation test of mutation-carriers vs wt
        re = st$cor$diff_test(cor_wt, cor_mut)
        list(mut = colnames(mut_hgnc),
             study = cor_row$tissue,
             path1 = pathways[1],
             path2 = pathways[2],
             effect = re$delta_cor[1,2],
             size = length(is_mut),
             p.value=re$p.value[1,2])
    }

    cur_mut = cur_mut[,colSums(cur_mut, na.rm=TRUE) >= 10, drop=FALSE]
    if (ncol(cur_mut) == 0)
        return(NULL)

    result = cur_mut %>%
        ar$split(along=2) %>%
        lapply(mut2assoc) %>%
        bind_rows() %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))
}

if (is.null(module_name())) {
    cor = io$load(INFILE) %>%
        mutate(adj.p = p.adjust(p.value,method="fdr")) %>%
        filter(adj.p < 0.01) %>%
        mutate(path1 = as.character(path1),
               path2 = as.character(path2),
               tissue = as.character(tissue))

    scores = io$load('../../scores/tcga/pathways_mapped/speed_matrix.RData')
    mut = io$load('../tcga_pathway_per_mutation/mut_driver_matrix.RData')
    cna = io$load('../tcga_pathway_per_mutation/cna_driver_matrix.RData')
    bem = ar$stack(list(mut, cna), along=2)
    bem = bem[substr(rownames(bem), 14, 16) == "01A",]
    ar$intersect(bem, scores, along=1)

    result = lapply(b$list$transpose(cor), do_assoc) %>%
        bind_rows()

    save(result, file=OUTFILE)
}
