library(dplyr)
io = import('io')

assocs = io$load('../../model/model_matrix.RData')$assocs

speed = io$load('../../model/model_matrix.RData')$assocs %>%
    group_by(pathway) %>%
    top_n(100, -p.value) %>%
	summarize(avg_fdr = mean(adj.p)) %>%
	print()
