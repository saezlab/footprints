library(dplyr)
b = import('base')
io = import('io')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/speed/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

scores = io$load(INFILE)
index = scores$index %>%
    mutate(id = rownames(.)) %>%
    select(id, pathway, cells, accession, effect) %>%
    arrange(pathway)

scores = scores$scores[index$id,]
scores[index$effect == "inhibiting",] = -scores[index$effect == "inhibiting",]
scores = t(scores)

rownames(index) = index$id
index$id = NULL

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

pheatmap::pheatmap(scores, annotation=index, scale="column", cluster_cols=FALSE)
