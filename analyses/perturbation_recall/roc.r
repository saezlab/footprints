library(cowplot)
library(dplyr)
library(reshape2)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
st = import('stats')
config = import('../../config')

#' Return a data.frame with precision-recall curves
#'
#' @param fid  Method ID descriptor ('speed_matrix', 'gsea_go', etc.)
methods2roc = function(fids) {
	data = lapply(fids, function(fid)
        io$load(module_file("../../scores/speed",
                paste0(fid, ".RData"), mustWork = TRUE)))

    lapply(data, scores2roc) %>%
        setNames(fids) %>%
        df$add_name_col(col="method", bind=TRUE) %>%
        na.omit()
}

#' Convert a scores matrix to TPR/FPR
#'
#' @param mat  Scores matrix [experiments x pathways]
scores2roc = function(mat) {
	index = mat$index # bit duplication from report/util_1 (4 lines w/ 2 above)
    sign = ifelse(index$effect == "activating", 1, -1)
	scores = mat$scores * sign # not sure if col+row-scale here, but we look per-pathway

    df = data.frame(pathway=index$pathway, scores, check.names=FALSE) %>%
        melt(id="pathway") %>%
        mutate(perturbed = as.integer(variable==pathway)) %>%
        select(-variable) %>%
        na.omit() %>%
        group_by(pathway) %>%
        do(st$roc(., "value", "perturbed")) %>%
        ungroup()
}

#' Create a ROC curve plot from a roc object
roc2plot = function(roc, width=1) {
    random_line = data.frame(x=c(0,1), y=c(0,1), method=roc$method[1])

    ggplot(roc, aes(x=FPR, y=TPR, color=method)) +
        geom_line(aes(x=x, y=y), data=random_line, color="grey", linetype="dashed", size=width) +
        geom_step(size=width) +
        coord_fixed() +
        facet_wrap(~pathway) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Calculate AUCs from a roc object
roc2auc = function(roc) {
    auc = roc %>%
        mutate(method = config$id2short(method)) %>%
        group_by(method, pathway) %>%
        arrange(value) %>%
        summarize(auc = st$roc_auc(value, perturbed==1)) %>%
        tidyr::spread(method, auc)
}

if (is.null(module_name())) {
    roc = methods2roc(setdiff(config$methods$analysis_set, "paradigm"))
    auc = roc2auc(roc)

    pdf("roc.pdf", paper="a4r", width=11, height=8)

    print(roc2plot(roc))

    dev.off()
}
