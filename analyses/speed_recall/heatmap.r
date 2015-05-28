library(dplyr)
b = import('base')
io = import('io')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/speed/speed_linear.r"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

scores = io$load(INFILE)
index = scores$index %>%
    mutate(id = rownames(.)) %>%
    select(id, pathway, cells, series) %>%
    arrange(pathway)
scores = scores$scores[index$id,]

pheatmap::pheatmap(scores, annotation=index, scale="column", cluster_cols=FALSE)
