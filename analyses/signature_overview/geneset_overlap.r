library(dplyr)
library(magrittr)
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
        select(gene, pathway) %>%
        group_by(pathway) %>% # unstack(gene ~ pathway) messes up JAK-STAT name
        tidyr::nest() %$%
        setNames(data, pathway) %>%
        lapply(function(x) x$gene) # nest() returns tibble otherwise

    names(sets) = paste0("gsva_", names(sets))
    sets = c(speed_matrix=list(speed), sets)
    names(sets) = config$id2short(names(sets))

    sets = lapply(seq_along(sets), function(i)
              data.frame(method = names(sets)[i], stack(sets[[i]]))) %>%
        bind_rows() %>%
        suppressWarnings() %>% # character conversion
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
    sorted = unique(sets$method)

    overlaps = sets %>%
        group_by(pathway) %>%
        do(set = set2overlap(.)) %>%
        ungroup() %>%
        tidyr::unnest() %>%
        mutate(method1 = factor(method1, levels=rev(sorted)),
               method2 = factor(method2, levels=sorted)) %>%
        filter(as.integer(method1) > length(sorted)-as.integer(method2))

    plt$matrix(overlaps, overlap ~ method1 + method2, palette="Blues") +
        geom_text(aes(label=overlap), size=3) +
        coord_fixed() +
        theme(legend.position = "none",
              text = element_text(size=10),
              axis.text.x = element_text(size=8),
              axis.text.y = element_text(size=8)) +
        facet_wrap(~pathway) + 
        xlab("") +
        ylab("")
}

if (is.null(module_name())) {
    pdf("geneset_overlap.pdf")
    on.exit(dev.off)

    sets = get_genesets()
    geneset_overlap_matrix(sets)
}
