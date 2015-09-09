library(dplyr)
b = import('base/util') # all of base: could not find function .local FIXME:
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')

fnames = list.files("../../scores/lincs", "\\.RData", full.names=TRUE)
names(fnames) = b$grep("/([a-zA-Z0-9_]+).RData", fnames)
index = io$load("../../util/lincs_perturbation_qc/index.RData") %>%
    filter(pert_type %in% c("control", "expression")) %>%
    as.data.frame()
index = index[!duplicated(index$distil_id),] #FIXME:
rownames(index) = index$distil_id
scores = lapply(fnames, io$load)

scores = ar$intersect_list(scores, along=1)
scores2 = ar$stack(scores, along=3) %>%
    ar$split(along=2, drop=TRUE) %>%
    lapply(function(x) ar$map(x, along=1, scale))
scores = scores2

# make sure index and scores have same rows
print(table(index$pathway))
scores$index = index
scores = ar$intersect_list(scores, along=1)
index = scores$index
scores$index = NULL
print(table(index$pathway))
scores = scores[unique(index$pathway)]

path2assocs = function(path) {
    message(path)
    score = scores[[path]]
    pathway = index$pathway == path
    score = score * ifelse(index$sign == "-", -1, 1) #FIXME:
    re = st$lm(score ~ pathway)
    re$pathway = path
    re
}
result = lapply(names(scores), path2assocs) %>%
    lapply(function(x) { #FIXME: this should be named properly already DPLYR BUG
        colnames(x) = sub("result\\.", "", colnames(x))
        x
    }) %>%
    bind_rows() %>%
    na.omit() %>%
    filter(term == "pathwayTRUE") %>%
    select(score, pathway, estimate, p.value)


# plot them as matrix
pdf("rna_perturb.pdf", paper="a4r")
on.exit(dev.off)
result %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"),
           label = ifelse(adj.p < 0.05, "*", ""),
           label = ifelse(adj.p < 0.005, "**", label),
           label = ifelse(adj.p < 1e-5, "***", label)) %>%
#    filter(adj.p < 0.1) %>%
    plt$matrix(estimate ~ pathway + score)
