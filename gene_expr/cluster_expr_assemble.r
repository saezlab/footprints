io = import('io')
clusters = io$load_regex("cluster_([A-Z]+)\\.RData")
save(clusters, file="expr_cluster.RData")
