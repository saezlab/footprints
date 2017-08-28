library(dplyr)
io = import('io')
st = import('stats')
util = import('./roc_util')

roc_df = function() {
    zdf = io$load(module_file('sigs_zscores.RData')) %>%
        util$scores2df() %>%
        mutate(inferred = sub("\\..*$", "", signature),
               method = "zscore")

    gsvadf = io$load(module_file('sigs_gsva.RData')) %>%
        util$scores2df() %>%
        mutate(inferred = sub("\\..*$", "", signature),
               method = "gsva")

#    fids = c("speed_matrix", "speed_original", "gsva_speed_matrix", "gsva_speed1")
    fids = c("speed_matrix", "gsva_speed_matrix", "speed_original", "speed_webserver", "epsa")
    asetdf = util$analysis_set(fids)

    # no scaling here because scores2df/analysis_set() take care of it
#    scoredf = dplyr::bind_rows(zdf, gsvadf, asetdf)
    scoredf = dplyr::bind_rows(zdf, asetdf)

    # per signature, how well do we infer perturbed pathway?
    roc = scoredf %>%
        na.omit() %>%
        mutate(matched = perturbed == inferred) %>%
        group_by(method, inferred, signature) %>%
        do(st$roc(., "score", "matched")) %>%
        ungroup()
}

plot_roc = function(roc, area=c('zscore','gsva')) {
    width=1
    random_line = data.frame(x=c(0,1), y=c(0,1), method=roc$signature[1])

    ggplot(roc, aes(x=FPR, y=TPR, color=method)) +
        geom_line(aes(x=x, y=y), data=random_line, color="grey", linetype="dashed", size=width) +
        stat_summary(data=filter(roc, method %in% area),
                     aes(fill=method), geom="ribbon", alpha=0.2,
                     fun.ymin = min, fun.ymax = max, color=NA) +
    #    geom_step(size=width, stat="summary", fun.y=median) +
    #    geom_step(size=0.1, alpha=0.1) +
        geom_step(data=filter(roc, ! method %in% area), size=width) +
        coord_fixed() +
        facet_wrap(~inferred) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

if (is.null(module_name())) {
    roc = roc_df()
    p = plot_roc(roc)
    auc = util$roc2auc(roc)

    pdf("roc_speed.pdf", paper="a4r", width=11, height=8)
    print(p)
    grid::grid.newpage()
    gridExtra::grid.table(auc)#, show.rownames=FALSE)
    dev.off()
}
