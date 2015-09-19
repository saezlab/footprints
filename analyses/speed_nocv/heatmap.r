library(dplyr)
b = import('base')
io = import('io')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/speed_train/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

data = io$load(INFILE)
index = data$index %>%
    select(id, pathway, effect) %>%
    arrange(pathway, effect) %>%
    as.data.frame()

scores = data$scores[index$id,]
scores = t(scores)
rownames(scores) = substr(rownames(scores), 0, 40)

# remove id column, add names for pheatmap to understand
rownames(index) = index$id
index$id = NULL

pdf(OUTFILE, paper="a4r", width=26, height=20)

pheatmap::pheatmap(scores,
                   annotation = index,
#                  scale = "column",
                   cluster_cols = FALSE,
                   show_colnames = FALSE,
                   annotation_legend = TRUE)

dev.off()
