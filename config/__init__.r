.io = import('io')

files = list.files(module_file(), "\\.yaml$")

for (.file in files)
    assign(sub("\\.yaml", "", .file), .io$read_yaml(.file))

if (is.null(module_name())) {
    query = commandArgs(TRUE)[1]
    path = strsplit(query, "/")[[1]]

    var = get(path[1])
    for (p in path[2:length(path)])
        var = var[[p]]

    cat(paste0(paste(var, collapse="\n"), "\n"))
}
