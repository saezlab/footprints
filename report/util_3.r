io = import('io')
vp = import('../analyses/drug_assocs/plot')

volcano = function(fid) {
    fp = io$file_path('assocs_mapped', fid, ext=".RData")
    assocs = vp$load_fun(fp)$assocs.pan
    vp$plot_pancan(assocs, text.size=2)
}

#' Returns names of a logical vector where the elements are TRUE
nst = function(x) names(x[x])

#' prepare stratifications of MAPK pathway
mapk_mut_pathway = function(scores, mut, pathway="MAPK", genes=c("BRAF","KRAS","NRAS")) {
    path = scores[names(tissues), "MAPK"]
    has_mut = apply(mut[,genes,drop=FALSE], 1, function(x) Reduce(`|`, x))
    path_active = path > quantile(path)[4] # top 25%
    path_inactive = path < quantile(path)[2] # bottom 25%
    path_null = path > quantile(path)[2] & path < quantile(path)[4] # rest

    list(
        # stratification by mutation
        "MAPK_all" = nst(!is.na(path)),
        "MAPK_wt" = nst(!has_mut),
        "MAPK_mut" = nst(has_mut),

        # stratification by SPEED score
        "MAPK+" = nst(path_active),
        "MAPK-" = nst(path_inactive),
        "MAPK_null" = nst(path_null),

        # stratification by SPEED score, wild-type subset
        "MAPK_wt" = nst(!has_mut),
        "MAPK+_wt" = nst(!has_mut & path_active),
        "MAPK-_wt" = nst(!has_mut & path_inactive),

        # stratification by SPEED score, mutated subset
        "MAPK_mut" = nst(has_mut),
        "MAPK+_mut" = nst(has_mut & path_active),
        "MAPK-_mut" = nst(has_mut & path_inactive)
    )
}

p53_mut_pathway = function(scores, mut, pathway="p53", genes="TP53") {
    path = scores[names(tissues), "p53"]
    has_mut = apply(mut[,genes,drop=FALSE], 1, function(x) Reduce(`|`, x))
    path_active = path > quantile(path)[4] # top 25%
    path_inactive = path < quantile(path)[2] # bottom 25%
    path_null = path > quantile(path)[2] & path < quantile(path)[4] # rest

    list(
        # stratification by mutation
        "TP53_all" = nst(!is.na(path)),
        "TP53_wt" = nst(!has_mut),
        "TP53_mut" = nst(has_mut),

        # stratification by SPEED score
        "p53+" = nst(path_active),
        "p53-" = nst(path_inactive),
        "p53_null" = nst(path_null),

        # stratification by SPEED score, wild-type subset
        "p53_wt" = nst(!has_mut),
        "p53+_wt" = nst(!has_mut & path_active),
        "p53-_wt" = nst(!has_mut & path_inactive),

        # stratification by SPEED score, mutated subset
        "p53_mut" = nst(has_mut),
        "p53+_mut" = nst(has_mut & path_active),
        "p53-_mut" = nst(has_mut & path_inactive)
    )
}
