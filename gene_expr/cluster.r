tissue2clusters = function(expr) {
    library(dplyr)
    ar = import('array')
    st = import('stats')

    # nmf-cluster them together, get optimal number of clusters
    # couple of hours w/o hpc
    clust = st$nmf(expr, k=2:10, max_iter=5000, rep=10) %>%
        arrange(k, cluster)
    coph = select(clust, k, rho) %>% unique()
    n = coph$k[which(diff(coph$rho) > 0) + 1]

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
#    hpc = import('hpc')

    TISSUE = commandArgs(TRUE)[1] %or% "BRCA"
    OUTFILE = commandArgs(TRUE)[2] %or% "cluster_BRCA.RData"

    dset = io$load("./corrected_expr.RData")
#    expr = ar$split(t(dset$corrected), along=2, subsets=dset$covar)
    #FIXME: along=1 should throw error
    #FIXME: along=2 should do the right thing

#    expr = list()
#    for (t in unique(dset$covar))
#       expr[[t]] = dset$corrected[,dset$covar==t]
#    clusters = lapply(expr, tissue2clusters)

    expr = dset$corrected[,dset$covar==TISSUE]
    rm(dset); gc()
    clusters = tissue2clusters(expr)

#    clusters = hpc$Q(tissue2clusters, expr, memory=10240)
    save(clusters, file=OUTFILE)
}
