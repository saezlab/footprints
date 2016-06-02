library(dplyr)
io = import('io')

summary = function() {
    speed = io$load(module_file('../../model/model_matrix.RData'))$assocs %>%
        group_by(pathway) %>%
        top_n(100, -p.value) %>%
        summarize(avg_fdr = mean(adj.p)) %>%
        transmute(Pathway = pathway, `Average FDR` = avg_fdr)
}

if (is.null(module_name())) {
	print(summary())
}
