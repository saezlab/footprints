library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
fdr_genes = import('../../analyses/signature_similarity/fdr_de')
path_map = import('../../util/genesets/mapped/speed1')$speed1 %>% unlist()

#' Send a single query to the SPEED webserver
#'
#' @param query  A character vector query set
query_geneset = function(query) {
    url = "http://speed.sys-bio.net/cgi-bin/run_speed.py#results"

    fill = list(
        genes = paste(head(query, 1000), collapse="\n"),
        zscore_percent = 1,
        overlap = 20,
        discard_percent = 50,
        unique_check = 1,
        dload = 1,
        background_file = "",
        run = "Submit"
    )
    fill = c(fill, setNames(path_map, rep("incl_pathways", length(path_map))))

    result = httr::POST(url, body=fill, encode="form")

    fdr = lapply(path_map, function(qq) as.numeric(b$grep(
        sprintf("#%s[\t0-9.]+\t([0-9.]+)\n", Hmisc::escapeRegex(qq)),
        httr::content(result)))) %>% unlist()
}

if (is.null(module_name())) {
    EXPR = commandArgs(TRUE)[2] %or% "../../data/expr.RData"
    ZSCORES = commandArgs(TRUE)[2] %or% "../../data/zscores.RData"
    OUTFILE = commandArgs(TRUE)[3] %or% "speed_webserver.RData"

    # load zscores, model building function, and expression for each experiment
    zdata = io$load(ZSCORES)
    index = zdata$index
    expr = io$load(EXPR)
    stopifnot(colnames(zdata$zscores) == names(expr$expr))

    # get the query sets for each experiment
    query_sets = mapply(fdr_genes$exp2sig, expr=expr$expr, index=expr$records) %>%
        lapply(names)
    query_sets = query_sets[sapply(query_sets, length) != 0]

    scores = lapply(query_sets, query_geneset) %>% simplify2array() %>% t()
    scores = -log10(scores)
    scores[is.na(scores)] = 0

    index = index[match(rownames(scores), index$id),]
    save(scores, index, file=OUTFILE)
}
