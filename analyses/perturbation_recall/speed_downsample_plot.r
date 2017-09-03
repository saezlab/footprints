library(ggplot2)
b = import('base')
io = import('io')

do_plot = function(aucdf) {
    ggplot(aucdf, aes(x=n_sigs, y=auc)) +
        geom_point() +
        stat_smooth() +
        facet_wrap(~ pathway, scales="free")
}

if (is.null(module_file())) {
    INFILE = commandArgs(TRUE)[1] %or% "speed_downsample_auc.RData"
    OUTFILE = commandArgs(TRUE)[1] %or% "speed_downsample_plot.pdf"

    aucdf = io$load(INFILE)
    pdf(OUTFILE)
    do_plot(aucdf)
    dev.off()
}
