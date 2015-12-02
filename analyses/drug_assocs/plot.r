library(dplyr)
b = import('base')
io = import('io')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

# load data
data = io$load(INFILE)
assocs.pan = data$assocs.pan
data$assocs.pan = NULL

# save pdf w/ pan-cancer & tissue specific
pdf(OUTFILE, paper="a4r", width=26, height=20)

# volcano plot for pan-cancer
assocs.pan %>%
    mutate(label = paste(Ys, scores, sep=":")) %>%
    plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
    plt$volcano(base.size=0.2) %>%
    print()

## matrix plot for pan-cancer
#assocs.pan %>%
#    mutate(lp = -log(adj.p),
#           label = ifelse(adj.p < 1e-2, '*', ''),
#           estimate = ifelse(adj.p < 0.1, estimate, NA)) %>%
#    plt$cluster(lp ~ scores + Ys, size=c(Inf,20)) %>%
#    plt$matrix(estimate ~ scores + Ys)

# volcano plot for tissue subsets
tissue_plot = function(assocs.tissue, name) assocs.tissue %>%
    mutate(label = paste(subset, Yf, scores, sep=":")) %>%
    plt$color$p_effect(pvalue="adj.p",
                       effect="estimate",
                       dir=-1,
                       thresh=0.1) %>%
    plt$volcano(p=0.1) + ggtitle(name)

for (x in seq_along(data))
    print(tissue_plot(data[[x]], names(data)[x]))

for (tissue in unique(data$assocs.tissue_noexp$subset)) {
    message(tissue) # if one fails easier to debug
    p = data$assocs.tissue_noexp %>%
        filter(subset == tissue) %>%
        mutate(label = paste(Yf, scores, sep=":")) %>%
        plt$color$p_effect(pvalue="adj.p",
                           effect="estimate",
                           dir=-1,
                           thresh=0.2) %>%
        plt$volcano(p=0.2) + ggtitle(tissue)
    print(p)
}

dev.off()
