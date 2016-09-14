library(dplyr)
b = import('base')
io = import('io')
plt = import('plot')

#' File loading helper for association objects
load_fun = function(fid, type="mapped") {
    io$load(module_file(paste0("assocs_", type), paste0(fid, ".RData")))
}

#' Volcano plot for pan-cancer
#'
#' @param assocs     data.frame of associations
#' @param base.size  scaling factor for point size
#' @param p          p-value cutoff for significance
#' @return           ggplot2 volcano object
plot_pancan = function(assocs, base.size=0.2, p=0.05, ...) {
    assocs %>%
        mutate(label = paste(drug, scores, sep=":")) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", thresh=p, dir=-1) %>%
        plt$volcano(base.size=base.size, p=p, ...)
}

#' Matrix plot for pan-cancer
#'
#' @param assocs  data.frame of associations
#' @return        ggplot2 matrix object
plot_matrix = function(assocs) {
    assocs %>%
        mutate(lp = -log(adj.p),
               label = ifelse(adj.p < 1e-2, '*', ''),
               estimate = ifelse(adj.p < 0.1, estimate, NA)) %>%
        plt$cluster(lp ~ scores + drug, size=c(Inf,20)) %>%
        plt$matrix(estimate ~ scores + drug)
}

#' Tissue-specific volcano plot
#'
#' @param assocs  data.frame of associations
#' @param name    title of the plot
#' @return        ggplot2 volcano object
plot_tissue = function(assocs, name) {
    assocs %>%
        mutate(label = paste(tissue, drug, scores, sep=":")) %>%
        plt$color$p_effect(pvalue = "adj.p",
                           effect = "estimate",
                           dir = -1,
                           thresh = 0.2) %>%
        plt$volcano(p = 0.2, label_top = 50) + ggtitle(name)
}

if (is.null(module_name())) {
    INFILE = commandArgs(TRUE)[1] %or% "assocs_mapped/speed_matrix.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "assocs_mapped/speed_matrix.pdf"

    # load data
    data = io$load(INFILE)

    # save pdf w/ pan-cancer & tissue specific
    pdf(OUTFILE, paper="a4r", width=11, height=8)

    print(plot_pancan(data$pan))
    print(plot_tissue(data$tissue$clinical, "clinical (min stage 2)"))
    print(plot_tissue(data$tissue$noexp, "no experimental (min stage 1, top 10 tissues)"))
    print(plot_tissue(data$tissue$sensi, "sensitive (5 lines measured, top 10 tissues)"))

    dev.off()
}
