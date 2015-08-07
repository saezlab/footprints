library(dplyr)
b = import('base/util') # all of base: could not find function .local FIXME:
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')

fnames = list.files("../../scores/lincs", "\\.RData", full.names=TRUE)
names(fnames) = b$grep("/([a-zA-Z0-9_]+).RData", fnames)
index = io$load("../../util/lincs_perturbation_qc/index.RData") %>%
    filter(sign %in% c("ctl_vehicle","trt_lig","trt_cp"))
index = index[!duplicated(index$distil_id),] #FIXME:
rownames(index) = index$distil_id
index = setNames(index$pathway, index$distil_id)
scores = lapply(fnames, io$load)

scores = ar$intersect_list(scores, along=1)
scores2 = ar$stack(scores, along=3) %>%
    ar$split(along=2, drop=TRUE) %>%
    lapply(function(x) ar$map(x, along=1, scale))
scores = scores2

scores = scores[!names(scores) %in% c("H2O2")]

result = lapply(names(scores), function(path) {
    message(path)
    score = scores[[path]]
    ar$intersect(score, index, along=1)
    pathway = as.matrix(index == path)
    re = st$lm(score ~ pathway)
    re$pathway = path
    re
}) %>%
    lapply(function(x) { #FIXME: this should be named properly already DPLYR BUG
        colnames(x) = sub("result\\.", "", colnames(x))
        x
    }) %>%
    do.call(dplyr::bind_rows, .) %>%
    na.omit() %>%
    filter(term == "pathwayTRUE") %>%
    select(score, pathway, estimate, p.value)


# plot them as matrix
pdf("matrix_plot.pdf", paper="a4r")
on.exit(dev.off)
result %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"),
           label = ifelse(adj.p < 0.05, "*", ""),
           label = ifelse(adj.p < 0.005, "**", label),
           label = ifelse(adj.p < 1e-5, "***", label)) %>%
#    filter(adj.p < 0.1) %>%
    plt$matrix(estimate ~ pathway + score)
