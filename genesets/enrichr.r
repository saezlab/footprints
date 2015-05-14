ll = import('base/list')
ar = import('array')

.enrichr_dir = file.path(module_file(), 'Enrichr/src/main/resources')

get_categories = function(enrichr_dir=.enrichr_dir) {
    doc = xmlInternalTreeParse(file.path(enrichr_dir, '/gene_set_libraries.xml'))
    categories = unname(xpathSApply(doc, "//categories/category/@name"))

    setNames(lapply(categories, function(f) unname(xpathSApply(doc, 
        paste0("//categories/category[@name='", f, "']/library/@name")))),
        make.names(categories))
}

parse_gmt = function(fname, weights=FALSE) {
    line2list = function(l) {
        re = strsplit(l, "\t")[[1]]
        targets = sapply(re[3:length(re)], function(f) strsplit(f, ","))
        genes = unname(sapply(targets, function(f) f[1]))
        if (weights) {
            if (length(targets[[1]]) > 1)
                val = sapply(targets, function(f) as.numeric(f[2]))
            else
                val = rep(1, length(targets))
            genes = setNames(val, genes)
        }
        setNames(list(toupper(genes)), re[1])
    }
    ll = readLines(fname)
    do.call(c, lapply(ll, line2list))
}

parse_category = function(lists, enrichr_dir=.enrichr_dir, weights=FALSE) {
    gmts = sapply(lists, function(l) paste0(file.path(enrichr_dir, l), ".gmt"))
    setNames(lapply(gmts, function(g) parse_gmt(g, weights=weights)), lists)
}
