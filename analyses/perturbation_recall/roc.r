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
#' @param scores  Scores object (scores: [experiments x pathways], index: df)
#' @param lookup  Lookup table for signature>pathway (default: colnames)
scores2roc = function(scores, lookup=setNames(colnames(scores$scores), colnames(scores$scores))) {
	index = scores$index # bit duplication from report/util_1 (4 lines w/ 2 above)
    sign = ifelse(index$effect == "activating", 1, -1)
	mat = scores$scores * sign # not sure if col+row-scale here, but we look per-pathway

    df = data.frame(perturbed = index$pathway,
                    mat,
                    check.names = FALSE) %>%
        tidyr::gather(signature, score, -perturbed) %>%
        mutate(inferred = lookup[signature],
               matched = as.integer(perturbed == inferred)) %>%
        na.omit() %>%
#        group_by(signature, inferred) %>%
#        group_by(perturbed) %>% #  current results
        group_by(signature) %>%
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
        facet_wrap(~signature) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Calculate AUCs from a roc object
roc2auc = function(roc) {
    auc = roc %>%
        mutate(method = config$id2short(method)) %>%
        group_by(method, signature) %>%
        arrange(score) %>%
        summarize(auc = st$roc_auc(score, matched==1)) %>%
        tidyr::spread(method, auc)
}

if (is.null(module_name())) {
    roc = methods2roc(setdiff(config$methods$analysis_set, "paradigm"))
#    auc = roc2auc(roc)

    pdf("roc.pdf", paper="a4r", width=11, height=8)

    print(roc2plot(roc))

    dev.off()
}
