io = import('io')
clusters = io$load_regex("\\.cluster_([A-Z]+)\\.RData", all.files=TRUE)
save(clusters, file="expr_cluster.RData")
