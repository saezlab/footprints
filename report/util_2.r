library(magrittr)
library(cowplot)
b = import('base')
io = import('io')
df = import('data_frame')
plt = import('plot')
config = import('../config')

#' Function to load data from the TCGA mutation associations
#'
#' @param dir  The subdir in the mutation folder
#' @param id   Identifier of the method
#' @return     data.frame with the association results
load_fun = function(dir, id) {
    fn = function(dir, id) {
        fname = paste0(id, ".RData")
        fpath = file.path("../analyses/tcga_pathway_per_mutation", dir, fname)
        io$load(module_file(fpath)) %catch% stop("File not found: ", fpath)
    }

    b$lnapply(id, function(id) fn(dir, id)) %>%
        df$add_name_col("method", bind=TRUE)
}

mut_assocs = load_fun("assocs_driver_mapped", config$methods$analysis_set)
mut_cov = mut_assocs %>%
    filter(subset == "pan_cov") %>%
    select(-subset)
mut_nocov = mut_assocs %>%
    filter(subset == "pan") %>%
    select(-subset)
cna_assocs = load_fun("assocs_cna_mapped", config$methods$analysis_set)
cna_cov = cna_assocs %>%
    filter(subset == "pan_cov") %>%
    select(-subset)
cna_nocov = cna_assocs %>%
    filter(subset == "pan") %>%
    select(-subset)

#' Plots a matrix with mutation-pathway associations for different methods
#'
#' @param assoc_obj  Data.frame from a linear model (st$lm)
#' @param gene       HGNC symbol of the gene (eg. TP53)
#' @return           ggplot object
mutation_method_plot = function(mut, assocs, p0=0.2, p1=0.05, p2=1e-3, ts=4) {
    p1 = filter(assocs, m == mut) %>%
        mutate(Wald = statistic) %>%
#        mutate(statistic = ifelse(adj.p < p0, statistic, NA)) %>%
        mutate(label = ifelse(adj.p < p1, "·", "")) %>%
        mutate(label = ifelse(adj.p < p2, "*", label)) %>%
        mutate(scores = config$pathways(scores, rev=TRUE),
               method = config$id2short(method)) %>%
        plt$matrix(Wald ~ scores + method, text_size=ts,
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
        mutate(m = b$refactor(m, -mlogp),
               method = config$id2short(method)) %>%
        ggplot(aes(x=m, y=mlogp, fill=method)) +
            geom_bar(stat="identity") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            geom_hline(aes(yintercept=-log10(0.05) * n_methods), linetype="dotted") +
            xlab("Mutated gene") +
            ylab("- log (p-value)")
}

#' Plots a matrix for associations of mutation/CNA and different methods
#'
#' @param n1,n2,n3,n4  Names of the aberration
#' @return             A cowplot arranged grid
mut_grid = function(n1, n2, n3, n4, ...) {
    name2obj = function(n) { # ifelse fails here and returns list
        if (grepl("_amp|_del", n))
            cna_cov
        else
            mut_cov
    }
    plot_grid(
        mutation_method_plot(n1, name2obj(n1), 0.2, 0.05, 1e-3, 10),
        mutation_method_plot(n2, name2obj(n2), 0.2, 0.05, 1e-3, 10),
        mutation_method_plot(n3, name2obj(n3), 0.2, 0.05, 1e-3, 10),
        mutation_method_plot(n4, name2obj(n4), 0.2, 0.05, 1e-3, 10),
        labels = c(n1, n2, n3, n4), ...
    ) +
      theme(axis.text = element_text(size=7),
            axis.title = element_text(size=9, face="bold"))
}
