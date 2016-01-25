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
