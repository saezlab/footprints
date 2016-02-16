library(magrittr)
library(cowplot)
b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
df = import('data_frame')
plt = import('plot')
tcga = import('data/tcga')

MUTFILE = "../analyses/tcga_pathway_per_mutation/mutations_annotated_pathwayactivities_v3_mikeformat.txt"

#' Loads a specific TCGA scores file
#' supplies path and extension, and cuts rownames at 16 chars
#'
#' @param x  An identifier, like 'speed_matrix'
#' @return   A matrix with samples as rows and genes and columns
load_fun = function(x) {
    re = io$file_path("../scores/tcga", x, ext=".RData") %>%
        io$load()
    rownames(re) = substr(rownames(re), 1, 16)
    re
}

#' Calculates pathway score associations per mutation
#'
#' @param s  A score matrix, sample IDs x pathways
#' @param m  A mutation matrix, sample IDs x genes
#' @return   A data.frame with the results for all pathways and mutations
lm_fun = function(s, m, study_covar=TRUE) {
    ar$intersect(m, s, along=1)
    study = tcga$barcode2study(rownames(m))
    if (study_covar)
        lm_result = st$lm(s ~ study + m)
    else
        lm_result = st$lm(s ~ m)

    lm_result %>%
        filter(!grepl("study", term)) %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))
}

#' Plots a matrix with mutation-pathway associations for different methods
#'
#' @param assoc_obj  Data.frame from a linear model (st$lm)
#' @param gene       HGNC symbol of the gene (eg. TP53)
#' @return           ggplot object
mutation_method_plot = function(assoc_obj, gene) {
    cur = filter(assoc_obj, m == gene)
    lmax = max(abs(cur$statistic))
    cur %>%
        mutate(statistic = ifelse(adj.p < 0.2, statistic, NA)) %>%
        mutate(label = ifelse(adj.p < 0.05, "*", "")) %>%
        mutate(label = ifelse(adj.p < 1e-3, "**", label)) %>%
        mutate(label = ifelse(adj.p < 1e-10, "***", label)) %>%
        plt$matrix(statistic ~ s + method, color="statistic", limits=c(-lmax,lmax)) +
        xlab("") + ylab("") + ggtitle(gene)
}

#' Plots stacked bar plots of overall association strength across different methods
#'
#' @param assoc_obj  Data.frame from a linear model (st$lm)
#' @param genes      HGNC symbols of the gene (character vector)
#' @return           ggplot object
mutation_overview_plot = function(assoc_obj, genes) {
    n_methods = length(unique(assoc_obj$method))

    assoc_obj %>%
        group_by(m, method) %>%
        summarize(mlogp = -sum(log10(adj.p))) %>%
        filter(m %in% genes) %>%
        ungroup() %>%
        mutate(m = b$refactor(m, -mlogp)) %>%
        ggplot(aes(x=m, y=mlogp, fill=method)) +
            geom_bar(stat="identity") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            geom_hline(aes(yintercept=-log10(0.05) * n_methods), linetype="dotted")
}

#' Plots matrices with mutation-pathway associations for different methods
#'
#' @param assoc_obj  Data.frame from a linear model (st$lm)
#' @param genes      HGNC symbols of the gene (character vector)
#' @param ...        Arguments passed to cowplot::plot_grid()
mutation_method_plots = function(assoc_obj, genes, ...) {
    plots = lapply(genes, function(x) mutation_method_plot(assoc_obj, x))
#    for (i in seq(1, length(plots), ncol)) { # iterate every 2 for page breaks
#        print(plot_grid(plotlist=plots[i:i+ncol], ncol=ncol, scale=scale))
#        cat("\n")
#    }
    print(plot_grid(plotlist=plots, ...))
}
