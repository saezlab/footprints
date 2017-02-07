library(dplyr)
library(magrittr)
b = import('base')
io = import('io')
config = import('../../config')

#' Returns gene sets of the different methods
#'
#' @param sets    Which gene sets to include (always includes speed_matrix)
#' @param mapped  Whether to use the pathway-mapped version (default: TRUE)
#' @param to_data_frame  Whether to convert results to a data_frame (default: TRUE)
#' @return        A data_frame of list of lists of different gene sets
get_genesets = function(sets = c("go", "reactome", "speed2016", "gatza"),
                        mapped = TRUE, to_data_frame=TRUE) {
    # import all sets from methods listed in config
    # missing: SPIA (KEGG gene sets, include?) and PARADIGM (no sets per se)
    if (mapped)
        path = "mapped"
    else
        path = "all"

    set_path = module_file(path)
    sets = b$lnapply(sets, function(s)
        io$load(file.path(set_path, paste0(s, ".RData"))))

    speed = io$load(module_file('../../model/model_matrix.RData'))$assocs %>%
        group_by(pathway) %>%
        top_n(100, -p.value) %>%
        ungroup() %>%
        select(gene, pathway) %>%
        group_by(pathway) %>% # unstack(gene ~ pathway) messes up JAK-STAT name
        tidyr::nest() %$%
        setNames(data, pathway) %>%
        lapply(function(x) x$gene) # nest() returns tibble otherwise

    names(sets) = paste0("gsva_", names(sets))
    sets = c(speed_matrix=list(speed), sets)
    names(sets) = config$id2short(names(sets))

    if (to_data_frame) {
        lapply(seq_along(sets), function(i)
               data.frame(method = names(sets)[i], stack(sets[[i]]))) %>%
            bind_rows() %>%
            suppressWarnings() %>% # character conversion
            select(method=method, pathway=ind, gene=values)
    } else
        sets
}
