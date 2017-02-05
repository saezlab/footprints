library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
plt = import('plot')
config = import('../../config')

get_genesets = function() { #FIXME: JAK.STAT
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

set2overlap = function(sets, fun=function(x,y) nrow(intersect(x,y))) {
    uniquify_colnames = function(df) {
        colnames(df) = make.names(colnames(df), unique=TRUE)
        df
    }

    sets %>%
        group_by(method) %>%
        tidyr::nest() %>%
        tidyr::crossing(., ., by="method") %>%
        uniquify_colnames() %>%
        transmute(method1 = method,
                  method2 = method.1,
                  overlap = purrr::map2_int(data, data.1, fun))
}

geneset_overlap_matrix = function(sets) {
    overlaps = sets %>%
        group_by(pathway) %>%
        do(set = set2overlap(.)) %>%
        ungroup() %>%
        tidyr::unnest()

    plt$matrix(overlaps, overlap ~ method1 + method2, palette="Blues") +
        geom_text(aes(label=all.vars(formula)[1])) +
        coord_fixed() +
        theme(legend.position = "none") +
        facet_wrap(~pathway) + 
        xlab("") +
        ylab("")
}

if (is.null(module_name())) {
    pdf("geneset_overlap.pdf", width=19, height=15)
    on.exit(dev.off)

    sets = get_genesets()
    geneset_overlap_matrix(sets)
}
