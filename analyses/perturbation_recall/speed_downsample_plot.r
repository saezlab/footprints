library(ggplot2)
b = import('base')
io = import('io')

aucdf = function(fname="speed_downsample_auc.RData") {
    io$load(module_file(fname))
}

do_plot = function(aucdf) {
    ggplot(aucdf, aes(x=n_sigs, y=auc)) +
        geom_point() +
        stat_smooth() +
        facet_wrap(~ pathway, scales="free")
}

if (is.null(module_name())) {
    INFILE = commandArgs(TRUE)[1] %or% "speed_downsample_auc.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "speed_downsample_plot.pdf"

    aucdf = aucdf(INFILE)
    pdf(OUTFILE)
    print(do_plot(aucdf))
    dev.off()
}
