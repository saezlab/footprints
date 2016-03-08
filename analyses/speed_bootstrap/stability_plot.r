do_plot = function() {
    library(dplyr)
    library(cowplot)
    io = import('io')
    scores = io$load('score_resample.RData')
        
    center = function(x) x - mean(x)
    var_mat = function(x) var(c(x))

    var_clines_given_bootstraps = apply(scores, c(1,2), center) %>%
        apply(3, var_mat)

    var_bootstraps_given_clines = apply(scores, c(2,3), center) %>%
        apply(2, var_mat)

    stability = var_bootstraps_given_clines / var_clines_given_bootstraps

    sdf = data.frame(pathway = names(stability), score = stability)

    ggplot(sdf, aes(x=pathway, y=score)) +
        geom_bar(stat="identity") +
        coord_flip() +
        scale_y_log10(breaks=c(1,2,5,10,100)) +
        geom_hline(yintercept=2, linetype="dashed") +
        geom_hline(yintercept=5, linetype="dashed") +
        xlab("Pathway") +
        ylab("Cell line variance over input variance")
}

if (is.null(module_name())) {
    pdf("stability.pdf")

    do_plot()

    dev.off()
}
