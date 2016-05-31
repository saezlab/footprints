library(dplyr)
b = import('base')
io = import('io')
plt = import('plot')

load_fun = function(fid) {
    io$load(module_file(fid))
}

plot_pancan = function(assocs.pan, base.size=0.2, p=0.05, ...) {
    # volcano plot for pan-cancer
    assocs.pan %>%
        mutate(label = paste(Ys, score, sep=":")) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", thresh=p, dir=-1) %>%
        plt$volcano(base.size=base.size, p=p, ...)
}

plot_matrix = function(assocs.pan) {
    # matrix plot for pan-cancer
    assocs.pan %>%
        mutate(lp = -log(adj.p),
               label = ifelse(adj.p < 1e-2, '*', ''),
               estimate = ifelse(adj.p < 0.1, estimate, NA)) %>%
        plt$cluster(lp ~ score + Ys, size=c(Inf,20)) %>%
        plt$matrix(estimate ~ score + Ys)
}

plot_tissue = function(assocs.tissue, name) {
    # volcano plot for tissue subsets
    assocs.tissue %>%
        filter(subset == name) %>%
        mutate(label = paste(tissue, drug, score, sep=":")) %>%
        plt$color$p_effect(pvalue="adj.p",
                           effect="estimate",
                           dir=-1,
                           thresh=0.1) %>%
        plt$volcano(p=0.1) + ggtitle(name)
}

if (is.null(module_name())) {
    INFILE = commandArgs(TRUE)[1] %or% "speed_matrix.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "speed_matrix.pdf"

    # load data
    data = io$load(INFILE)

    # save pdf w/ pan-cancer & tissue specific
    pdf(OUTFILE, paper="a4r", width=11, height=8)

    print(plot_pancan(data$assocs.pan))

    print(plot_tissue(filter(data, subset=="clinical"), "clinical"))
    print(plot_tissue(filter(data, subset=="noexp"), "no experimental"))
    print(plot_tissue(filter(data, subset=="sensi"), "sensitive"))

    dev.off()
}
