library(dplyr)
library(cowplot)
io = import('io')

load_fun = function(fname, name) {
    re = io$load(module_file(fname)) %>%
        select(id, pathway, x, y)
    re$method = name
    re
}

do_plot = function() {    
    df = bind_rows(load_fun('tsne_fc.RData', "fold change"),
                   load_fun('tsne_fullmat.RData', "no filter"),
                   load_fun('tsne_top100.RData', "top 100 genes"))

    ggplot(df, aes(x=x, y=y, color=pathway)) +
        geom_point() +
        facet_wrap(~method) +
        xlab("A.U.") +
        ylab("A.U.") +
        coord_fixed() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank())

}

if (is.null(module_name())) {
    OUTFILE = "plot.pdf"

    pdf(OUTFILE, paper="a4r", width=10, height=4)
    do_plot()
    dev.off()
}
