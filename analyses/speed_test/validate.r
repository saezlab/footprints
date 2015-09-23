library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
df = import('data_frame')

# load all in ../scores/speed/*.RData
objs = io$load_regex("(.*)\\.RData", "../../scores/speed_test")
index = objs$speed_linear$index # index is the same on all objs
# scale the pathway scores for each sample to mean=0, sd=1
sfun = function(o)
    ar$map(o$scores, along=2, function(x) as.numeric(x)) # scale here??
scores = lapply(objs, sfun) %>% ar$intersect_list(along=1)
index = index[match(rownames(scores[[2]]), index$id),]
#TODO: why dropping one here?

# make scores a list of pathways, with matrix cols = methods
test_hyper = function() {
    method2concordance = function(path, method) {
        sub_idx = filter(index, path == pathway)
        score = scores[[method]][sub_idx$id,]

        reverse = ifelse(sub_idx$effect == "activating", FALSE, TRUE)
        con_fun = function(s,r) names(sort(s,decreasing=!r))[1] == path
        concordance = mapply(con_fun, ar$split(score, along=1, drop=TRUE), reverse) %>%
            b$omit$na()

        1.0 - phyper(q = sum(concordance) - 1,
                     m = length(concordance),
                     n = length(concordance) * ncol(score),
                     k = length(concordance))
    }

    idx_df = df$create_index(path = unique(index$pathway),
                             method = names(scores),
                             expand_grid = TRUE)
    df$call(idx_df, method2concordance) %>%
        mutate(adj.p = p.adjust(result, method="fdr"),
               logp = log(result))
}

result = ar$construct(result ~ path + method, data=test_hyper())
pdf("validate.pdf", paper="a4r", width=26, height=20)
pheatmap::pheatmap(result,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE)
dev.off()
