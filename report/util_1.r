library(dplyr)
library(gridExtra)
library(corrplot)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
config = import('../config')

#' Takes a method subset and plots the pathway associations and individual scores
#'
#' @param fid    A file identifer, e.g. 'speed_matrix' or 'gsea_reactome'
#' @param title  Title of the plot
#' @return       A ggplot matrix plot and heatmap in one line
perturb_score_plots = function(fid, title=NULL) {
    data = io$file_path("../scores/speed", fid, ext=".RData") %>%
        io$load()

	index = data$index
	scores = data$scores
	sign = ifelse(index$effect == "activating", 1, -1)
	pathway = sign * ar$mask(index$pathway) + 0

	# the associations per pathway
	result = st$lm(scores ~ pathway) %>%
		mutate(mlogp = -log(p.value)) %>%
		mutate(label = ifelse(p.value < 1e-5, "·", "")) %>%
		mutate(label = ifelse(p.value < 1e-10, "*", label)) %>%
        mutate(pathway = config$pathways(pathway),
               scores = config$pathways(scores, rev=TRUE))

	p1 = plt$matrix(result, statistic ~ scores + pathway, palette="RdBu", symmetric=TRUE, text_size=6) +
        theme(axis.text = element_text(size = 8)) +
		xlab("Pathway perturbed") +
		ylab("Assigned score")
        coord_fixed()

	# and individual experiments
	annot = data$index %>%
        arrange(pathway) %>%
		select(id, pathway, effect) %>%
		arrange(pathway, effect) %>%
		as.data.frame()

	scores = data$scores[annot$id,]
	scores[annot$effect == "inhibiting"] = - scores[annot$effect == "inhibiting"]
	scores = t(scores)
	rownames(scores) = substr(rownames(scores), 0, 40)

	# remove id column, add names for pheatmap to understand
	rownames(annot) = annot$id
	annot$id = NULL

	p2 = pheatmap::pheatmap(scores,
                            annotation = annot,
#                            scale = "column",
                            cluster_cols = FALSE,
                            show_colnames = FALSE,
                            annotation_legend = FALSE,
                            treeheight_row = 20,
                            legend = FALSE,
                            cellwidth = 0.3, silent=TRUE)

    grid.arrange(p1, p2$gtable, ncol=2, top=title)
}

cor_plots = function(fid) {
    gdsc = io$file_path("../scores/gdsc/pathways_mapped", fid, ext=".RData") %>%
        io$load()

    tcga = io$file_path("../scores/tcga/pathways_mapped", fid, ext=".RData") %>%
        io$load() %>%
        tcga$filter(primary=TRUE, cancer=TRUE) %>%
        na.omit()

    par(mfrow=c(1, 2))
    corrplot(cor(tcga), title="TCGA primary tumors")
    corrplot(cor(gdsc), title="GDSC cell lines")
    par(mfrow=c(1, 1))
}
