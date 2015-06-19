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

    INFILE = commandArgs(TRUE)[1] %or% "../model/model_linear.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "speed_cluster.RData"

    dset = io$load("./corrected_expr.RData")
    zfit = io$load(INFILE)
    expr = dset$corrected

    ar$intersect(zfit, expr, along=1)
    scores = t(expr) %*% zfit %>%
        ar$map(along=1, base::scale) # pathway across experiments

#    expr = ar$split(t(dset$corrected), along=2, subsets=dset$covar)
    #FIXME: along=1 should throw error
    #FIXME: along=2 should do the right thing
    expr = list()
    for (t in unique(dset$covar))
       expr[[t]] = t(scores[dset$covar==t,])
    clusters = lapply(expr, tissue2clusters)

    save(clusters, file=OUTFILE)
}
