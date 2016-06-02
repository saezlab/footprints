.b = import('base')
.io = import('io')

files = list.files(module_file(), "\\.yaml$", full.names=TRUE) %>%
    setNames(.b$grep("/([^/]+)\\.yaml", .))

for (.i in seq_along(files))
    assign(names(files)[.i], .io$read_yaml(files[.i]))

if (is.null(module_name())) {
    query = commandArgs(TRUE)[1]
    path = strsplit(query, "/")[[1]]

    var = get(path[1])
    for (p in path[2:length(path)])
        var = var[[p]]

    cat(paste0(paste(var, collapse="\n"), "\n"))
}

alph = function(x) rev(gtools::mixedsort(unique(x)))

id2name = function(ids, drop=TRUE) {
    idx = match(ids, methods$ids)
    fac = factor(methods$names[idx], levels=methods$names)
    if (drop)
        droplevels(fac)
    else
        factor
}

id2short = function(ids, drop=TRUE) {
    idx = match(ids, methods$ids)
    fac = factor(methods$short[idx], levels=methods$short)
    if (drop)
        droplevels(fac)
    else
        factor
}
