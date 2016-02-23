library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/speed/gsea_reactome.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

data = io$load(INFILE)
index = data$index %>%
    select(id, pathway, effect) %>%
    arrange(pathway, effect) %>%
    as.data.frame()

# scale each pathway across samples, then pathways per sample
scores = data$scores[index$id,]
scores[index$effect == "inhibiting"] = - scores[index$effect == "inhibiting"]
scores = t(scores)
rownames(scores) = substr(rownames(scores), 0, 40)

# remove id column, add names for pheatmap to understand
rownames(index) = index$id
index$id = NULL

# this does the same as scale columns
scale_pheatmap = function(mat, n=100, center=TRUE) {
    color = 1:n
#    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(n)

    # scale matrix
    x = ar$map(mat, along=1, scale)

    # generate breaks
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    breaks = seq(-m, m, length.out = n + 1)

    scaled = color[as.numeric(cut(c(x), breaks=breaks, include.lowest=TRUE))]
    mat[] = scaled
    mat
}

pdf(OUTFILE, paper="a4r", width=26, height=20)

pheatmap::pheatmap(scores,
                   annotation = index,
                   cluster_cols = FALSE,
                   show_colnames = FALSE,
                   annotation_legend = TRUE)

# this does the same as scale columns
#pheatmap::pheatmap(scale_pheatmap(scores),
#                   annotation = index,
#                   cluster_cols = FALSE,
#                   show_colnames = FALSE,
#                   annotation_legend = TRUE)

pheatmap::pheatmap(scores,
                   annotation = index,
                   scale = "column",
                   cluster_cols = FALSE,
                   show_colnames = FALSE,
                   annotation_legend = TRUE)

dev.off()
