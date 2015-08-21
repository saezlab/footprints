io = import('io')
ar = import('array')
st = import('stats')
tcga = import('data/tcga')
plt = import('plot')

mut = io$load('mut_filtered.RData')

fnames = list(
    speed = "../../scores/tcga/speed_linear.RData",
    reactome = "../../scores/tcga/reactome.RData",
    go = "../../scores/tcga/go.RData",
    pathifier = "../../scores/merge/pathifier.RData",
    spia = "../../scores/merge/spia.RData"
)
scores = lapply(fnames, function(fname) {
    re = io$load(fname)
    rownames(re) = substr(rownames(re), 1, 16)
    re
})

scores = ar$intersect_list(scores, along=1)
colnames(scores$spia) = c("Trail","NFkB","EGFR","TGFb","VEGF","MAPK","PI3K","JAK-STAT") #FIXME:
scores2 = ar$stack(scores, along=3) %>%
    ar$split(along=2, drop=TRUE) %>%
    lapply(function(x) ar$map(x, along=1, scale))
scores = scores2[names(scores2) != "H2O2"]

mut = filter(mut, barcode %in% rownames(scores[[1]]))
mut$path = as.character(mut$path) # fix stringasfactors!
mutations = matrix(FALSE, nrow=nrow(scores[[1]]), ncol=length(scores),
                   dimnames=list(rownames(scores[[1]]), names(scores)))
for (path in names(scores))
    mutations[mut$barcode[mut$path==path], path] = TRUE

tissues = as.matrix(tcga$barcode2index(rownames(mutations))$Study.Abbreviation)
#TODO: this should work
#result = st$lm(scores ~ tissues + mutations, group=c("pathways","mutations"))
#    filter(term == "mutationsTRUE") %>%
#    select(-term, -tissues, -size)

result = lapply(names(scores), function(path) {
    score = scores[[path]]
    mutation = mutations[,path,drop=FALSE]
    st$lm(score ~ tissues + mutation)
}) %>%
    lapply(function(x) { #FIXME: this should be named properly already
        colnames(x) = sub("result\\.", "", colnames(x))
        x
    }) %>%
    do.call(dplyr::bind_rows, .) %>%
    na.omit() %>%
    filter(term == "mutationTRUE") %>%
    select(score, mutation, estimate, p.value)


# plot them as matrix
pdf("matrix_plot.pdf", paper="a4r")
on.exit(dev.off)
result %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"),
           label = ifelse(adj.p < 0.05, "*", ""),
           label = ifelse(adj.p < 0.005, "**", label),
           label = ifelse(adj.p < 1e-5, "***", label)) %>%
#    filter(adj.p < 0.1) %>%
    plt$matrix(estimate ~ mutation + score)
