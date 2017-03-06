library(cowplot)
library(dplyr)
library(reshape2)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
st = import('stats')
config = import('../../config')

#' Creates a data.frame for a scores object
#'
#' @param sobj  A file name or scores object
scores2df = function(sobj) {
    if (is.character(sobj))
        sobj = io$load(sobj)

    sign = ifelse(sobj$index$effect == "activating", 1, -1)
    scores = sobj$scores * sign
#    scores = ar$map(scores, along=2, function(x) scale(x,center=F))

    df = data.frame(perturbed = sobj$index$pathway,
                    scores,
                    check.names = FALSE) %>%
        tidyr::gather(signature, score, -perturbed)
}

#' Returns the analysis set and their relative paths
analysis_set2roc = function() {
    fids = setdiff(config$methods$analysis_set, "paradigm")

    # data.frame with columns:
    #   perturbed : which pathway was perturbed in the experiment
    #   signature : which signature was used to quantify
    #   score     : the score assigned by the signature
    #   method    : which method the signature comes from
    scoredf = paste0(fids, ".RData") %>%
        module_file("../../scores/speed", ., mustWork = TRUE) %>%
        setNames(fids) %>%
        lapply(scores2df) %>%
        df$add_name_col(col="method", bind=TRUE)

    # group_by:
    #   perturbed : are scores higher for pathway sig than others?
    #   signature : does sig assign higher score for pathway that is perturbed?
    roc = scoredf %>%
        na.omit() %>%
        mutate(matched = perturbed == signature) %>%
        group_by(perturbed, method) %>%
        do(st$roc(., "score", "matched")) %>%
        ungroup()
}

#' Create a ROC curve plot from a roc object
roc2plot = function(roc, width=1) {
    random_line = data.frame(x=c(0,1), y=c(0,1), method=roc$method[1])

    ggplot(roc, aes(x=FPR, y=TPR, color=method)) +
        geom_line(aes(x=x, y=y), data=random_line, color="grey", linetype="dashed", size=width) +
        geom_step(size=width) +
        coord_fixed() +
        facet_wrap(~perturbed) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Calculate AUCs from a roc object
roc2auc = function(roc) {
    auc = roc %>%
        mutate(method = config$id2short(method)) %>%
        group_by(method, perturbed) %>%
        arrange(score) %>%
        summarize(auc = st$roc_auc(score, matched==1)) %>%
        tidyr::spread(method, auc)
}

if (is.null(module_name())) {
    roc = analysis_set2roc()
#    auc = roc2auc(roc)

    pdf("roc.pdf", paper="a4r", width=11, height=8)

    print(roc2plot(roc))

    dev.off()
}
