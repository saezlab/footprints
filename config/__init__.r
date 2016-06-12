library(ggplot2)
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

.rev = function(x, rev=TRUE) {
    if (rev)
        base::rev(x)
    else
        x
}

pathways = function(ids, rev=FALSE) {
    factor(ids, levels=.rev(gtools::mixedsort(unique(ids)), rev=rev))
}

id2name = function(ids, drop=TRUE, rev=FALSE) {
    idx = match(ids, c(methods$ids, methods$additional_name_mapping$ids))
    name = c(methods$names, methods$additional_name_mapping$names)
    fac = factor(name[idx], levels=.rev(name, rev=rev))
    if (drop)
        droplevels(fac)
    else
        factor
}

id2short = function(ids, drop=TRUE, rev=FALSE) {
    idx = match(ids, c(methods$ids, methods$additional_name_mapping$ids))
    short = c(methods$short, methods$additional_name_mapping$short)
    fac = factor(short[idx], levels=.rev(short, rev=rev))
    if (drop)
        droplevels(fac)
    else
        factor
}

facet_theme = theme(panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					strip.background = element_blank(),
					panel.border = element_blank())
