tissue2clusters = function(expr) {
    library(dplyr)
    ar = import('array')
    st = import('stats')

    # nmf-cluster them together, get optimal number of clusters
    # couple of hours w/o hpc
    clust = st$nmf(expr, k=2:10, max_iter=10000, rep=10) %>%
        arrange(k, cluster)
    coph = select(clust, k, rho) %>% unique()
    n = coph$k[which(diff(sign(diff(coph$rho)))==-2)+1]

    # create ar$mask for tcga and gdsc separately
    lapply(n, function(n) {
        df = clust %>% filter(k == n)
        re = setNames(df$cluster, df$sample) %>% ar$mask() + 0
        colnames(re) = paste0("k", n, "_", colnames(re))
        re
    })
}

if (is.null(module_name())) {
    b = import('base')
    io = import('io')
    ar = import('array')

    TISSUE = commandArgs(TRUE)[1] %or% "BRCA"
    OUTFILE = commandArgs(TRUE)[2] %or% "cluster_BRCA.RData"

    tissue = io$h5load("corrected_expr.h5", "/tissue")
    expr = t(io$h5load("corrected_expr.h5", "/expr", 
                       index = which(tissue==TISSUE)))

    clusters = tissue2clusters(expr)
    save(clusters, file=OUTFILE)
}
