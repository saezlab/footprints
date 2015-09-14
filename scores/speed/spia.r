library(dplyr)
io = import('io')
ar = import('array')

pathifier = io$load("../../data/spia_scores.RData")
scores = pathifier$scores[!is.na(pathifier$scores)]
records = pathifier$records[!is.na(pathifier$scores)]

names(scores) = 1:length(scores)
scores = ar$stack(scores, along=1)
records = records[as.integer(rownames(scores))]

index = lapply(records, function(r)
    r[c('accession', 'platform', 'pathway', 'cells', 'treatment', 'effect', 'hours')]
) %>% setNames(., 1:length(.)) %>% bind_rows() %>% as.data.frame()

save(index, scores, file="spia.RData")
