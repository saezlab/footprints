library(dplyr)
library(treemap)
io = import('io')

index = io$load('../../data/zscores.RData')$index

tmap = c(table(index$pathway))
tmap = data.frame(pathway = names(tmap), n = unname(tmap)) %>%
    mutate(pathway = sprintf("%s (%i)", pathway, n))

pdf('treemap.pdf')
treemap(tmap, 'pathway', 'n')
dev.off()
