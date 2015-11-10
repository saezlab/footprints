library(ggplot2)
library(reshape2)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
st = import('stats')
tcga = import('data/tcga')

PATHWAY = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_linear.RData"
COR = commandArgs(TRUE)[2] %or% "correlation.RData"
MUTFILE = "../tcga_pathway_per_mutation/mutations_annotated_pathwayactivities_v3_mikeformat.txt"
OUTFILE = commandArgs(TRUE)[2] %or% "assocs.pdf"

#' @param subs  The subset or study ("pan" or TCGA cancer type identifier)
#' @param 
do_assoc = function(cor_row) {
    cur_mut = mut
    if (cor_row$tissue != "pan")
        cur_mut = filter(cur_mut, study == cor_row$tissue)

    mut2assoc = function(mut_hgnc) {
#        message(mut_hgnc)
        score = paths[,c(cor_row$path1, cor_row$path2)]

        has_mut = filter(cur_mut, hgnc==mut_hgnc) %>%
            select(sample) %>%
            unlist(use.names=FALSE) %>%
            intersect(rownames(score))
        has_wt = mut$sample[! mut$sample %in% has_mut] %>%
            intersect(rownames(score))

        # do correlation test of mutation-carriers vs wt
        re = st$cor$diff_test(score[has_wt,], score[has_mut,], return_effect=TRUE)
        list(label = paste(c(mut_hgnc, colnames(score)), collapse=":"),
             effect = re$effect[1,2],
             p.value=re$p.value[1,2])
    }

    cur_mut %>%
        group_by(hgnc) %>%
        filter(n() >= 10) %>%
        ungroup() %>%
        select(hgnc) %>%
        unlist(use.names=FALSE) %>%
        unique() %>%
        lapply(function(x) mut2assoc(x) %catch% NULL) %>%
        bind_rows() %>%
        mutate(study = cor_row$tissue)
}

cor = io$load(COR) %>%
    mutate(adj.p = p.adjust(p.value,method="fdr")) %>%
    filter(adj.p < 0.01) %>%
    mutate(path1 = as.character(path1),
           path2 = as.character(path2),
           tissue = as.character(tissue))

mut = io$read_table(MUTFILE, header=TRUE) %>%
    transmute(hgnc = GENE_NAME,
              sample = substr(Tumor_Sample_Barcode, 1, 16),
              study = tcga$barcode2study(Tumor_Sample_Barcode)) %>%
    filter(!is.na(study) & study != "READ")

paths = io$load(PATHWAY)

re = lapply(b$list$transpose(cor), do_assoc) %>%
    bind_rows() %>%
    group_by(study) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))
