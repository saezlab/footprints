library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')

perturb_score_plots = function(fid) {
#	library(gridGraphics)
	library(grid)
    library(gridBase)
#	library(gridExtra)
    par(mfrow=c(1,2))

    data = io$file_path("../scores/speed", fid, ext=".RData") %>%
        io$load()

	index = data$index
	scores = data$scores
	sign = ifelse(index$effect == "activating", 1, -1)
	pathway = sign * ar$mask(index$pathway) + 0

	# the associations per pathway
	result = st$lm(scores ~ pathway) %>%
		mutate(mlogp = -log(p.value)) %>%
		mutate(label = ifelse(p.value < 1e-5, "*", "")) %>%
		mutate(label = ifelse(p.value < 1e-10, "***", label))

	p1 = plt$matrix(result, statistic ~ scores + pathway, palette="RdBu", symmetric=TRUE) +
		xlab("Pathway perturbed") +
		ylab("Assigned score")

    plot.new()
    vps = baseViewports()
    pushViewport(vps$figure)
    vp1 = plotViewport(c(0,0,0,0))
    print(p1, vp=vp1)
    popViewport()

	# and individual experiments
	annot = data$index %>%
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

    plot.new()
	pheatmap::pheatmap(scores,
					   annotation = annot,
					   scale = "column",
					   cluster_cols = FALSE,
					   show_colnames = FALSE,
					   annotation_legend = FALSE,
                       treeheight_row = 20,
                       legend = FALSE,
                       cellwidth = 0.3)
}
