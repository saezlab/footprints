library(magrittr)
library(cowplot)
b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
df = import('data_frame')
plt = import('plot')
tcga = import('data/tcga')

#' Loads a specific TCGA scores file
#'
#' supplies path and extension, and cuts rownames at 16 chars
#'
#' @param x       An identifier, like 'speed_matrix'
#' @param filter  Which subset of the mutation associations to use
#' @return        A matrix with samples as rows and genes and columns
assoc_df = function(x, filter="pan_cov") {
    load_fun = function(x) {
        io$file_path("../analyses/tcga_pathway_per_mutation/assocs_driver_mapped", x, ext=".RData") %>%
            io$load() %>%
            filter(subset == filter) %>%
            mutate(method = x)
    }
    re = bind_rows(lapply(x, load_fun))
}

scores = c('gsea_reactome', 'gsea_go', 'speed_matrix', 'pathifier', 'spia', 'paradigm')

assocs_cov = assoc_df(scores, "pan_cov")
assocs_nocov = assoc_df(scores, "pan")

alph = function(x) rev(gtools::mixedsort(unique(x)))

#' Plots a matrix with mutation-pathway associations for different methods
#'
#' @param assoc_obj  Data.frame from a linear model (st$lm)
#' @param gene       HGNC symbol of the gene (eg. TP53)
#' @return           ggplot object
mutation_method_plot = function(mut, assocs, p0=0.2, p1=0.05, p2=1e-3) {
    p1 = filter(assocs, m == mut) %>%
        mutate(Wald = statistic) %>%
#        mutate(statistic = ifelse(adj.p < p0, statistic, NA)) %>%
        mutate(label = ifelse(adj.p < p1, "·", "")) %>%
        mutate(label = ifelse(adj.p < p2, "*", label)) %>%
        mutate(scores = factor(scores, levels=alph(scores)),
               method = factor(method, levels=alph(method))) %>%
        plt$matrix(Wald ~ scores + method, text_size=4,
                   color="Wald", symmetric=TRUE, reverse_colors=TRUE) +
        xlab("") + ylab("")
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
            geom_hline(aes(yintercept=-log10(0.05) * n_methods), linetype="dotted") +
            xlab("Mutated gene") +
            ylab("- log (p-value)")
}
