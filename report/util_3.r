io = import('io')
vp = import('../analyses/drug_assocs/plot')

volcano = function(fid) {
    fp = io$file_path('assocs_mapped', fid, ext=".RData")
    assocs = vp$load_fun(fp)$assocs.pan
    vp$plot_pancan(assocs, text.size=2)
}

#' Returns names of a logical vector where the elements are TRUE
nst = function(x) names(x[x])

#' Generates a stratification structure using pathway and mutations
#'
#' @param pathway    Our pathway descriptor, e.g. "MAPK"
#' @param mutations  A character vector of which genes should be
#'                   used as markers to compare to
#' @param strat      For stratifications, look for subsets in "mut" or "wt"?
stratify = function(scores, mut, pathway="MAPK", genes=c("BRAF","KRAS","NRAS")) {
    path = scores[names(tissues), pathway]
    has_mut = apply(mut[, genes, drop=FALSE], 1, any)
    path_active = path > quantile(path)[4] # top 25%
    path_inactive = path < quantile(path)[2] # bottom 25%
    path_null = path > quantile(path)[2] & path < quantile(path)[4] # rest

    re = list(
        # stratification by mutation
        "GENE_all" = nst(!is.na(path)),
        "GENE_wt" = nst(!has_mut),
        "GENE_mut" = nst(has_mut),

        # stratification by SPEED score
        "PATH+" = nst(path_active),
        "PATH-" = nst(path_inactive),
        "PATH_null" = nst(path_null),

        # stratification by SPEED score, wild-type subset
        "PATH_wt" = nst(!has_mut),
        "PATH+_wt" = nst(!has_mut & path_active),
        "PATH-_wt" = nst(!has_mut & path_inactive),

        # stratification by SPEED score, mutated subset
        "PATH_mut" = nst(has_mut),
        "PATH+_mut" = nst(has_mut & path_active),
        "PATH-_mut" = nst(has_mut & path_inactive)
    )

    names(re) = sub("PATH", pathway, names(re))
    if (length(genes) == 1)
        names(re) = sub("GENE", genes, names(re))
    else
        names(re) = sub("GENE", pathway, names(re))
    re
}
