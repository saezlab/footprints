library(cowplot)
library(dplyr)
library(reshape2)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
#config = io$load('../../config')# FIXME:
config = c("speed_matrix", "gsea_go", "gsea_reactome", "spia", "pathifier")

#' Return a data.frame with precision-recall curves
#'
#' @param fid  Method ID descriptor ('speed_matrix', 'gsea_go', etc.)
method2pr_df = function(fid) {
	data = io$file_path("../../scores/speed", fid, ext=".RData") %>%
        io$load()

	index = data$index # bit duplication from report/util_1 (4 lines w/ 2 above)
    sign = ifelse(index$effect == "activating", 1, -1)
	scores = data$scores * sign # not sure if col+row-scale here, but we look per-pathway

    df = data.frame(pathway=index$pathway, scores) %>%
        melt(id="pathway") %>%
        mutate(perturbed = as.integer(variable==pathway),
               method = fid) %>%
        select(-variable) %>%
        group_by(pathway) %>%
        do(st$roc(., "value", "perturbed")) %>%
        ungroup()
}

roc = lapply(config, method2pr_df) %>%
    bind_rows()

pdf("roc.pdf", paper="a4r", width=11, height=8)

ggplot(roc, aes(x=FPR, y=TPR, color=method)) +
    geom_step() +
    facet_wrap(~pathway)

dev.off()
