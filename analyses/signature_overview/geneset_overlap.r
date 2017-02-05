library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
plt = import('plot')
config = import('../../config')

get_genesets = function() {
    # import all sets from methods listed in config
    # missing: SPIA (KEGG gene sets, include?) and PARADIGM (no sets per se)
    sets = c("go", "reactome", "speed2016", "gatza")

    set_path = module_file('../../util/genesets/mapped')
    sets = b$lnapply(sets, function(s)
        io$load(file.path(set_path, paste0(s, ".RData"))))

    speed = io$load(module_file('../../model/model_matrix.RData'))$assocs %>%
        group_by(pathway) %>%
        top_n(100, -p.value) %>%
        ungroup() %>%
        unstack(gene ~ pathway) %>%
        ar$split(along=2, drop=TRUE)

    names(sets) = paste0("gsva_", names(sets))
    sets = c(speed_matrix=list(speed), sets)
    names(sets) = config$id2short(names(sets))

    sets = lapply(seq_along(sets), function(i)
              data.frame(method = names(sets)[i], stack(sets[[i]]))) %>%
        bind_rows() %>%
        select(method=method, pathway=ind, gene=values)
}

geneset_overlap_matrix = function(sets) {
    # compute their overlap

    # return a ggplot upper diagonal matrix
}

if (is.null(module_name())) {
    pdf("geneset_overlap.pdf", width=19, height=15)
    on.exit(dev.off)

    sets = get_genesets()
    set_overlap_matrix(sets)
}
