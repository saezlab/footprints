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

test_wilcox = function() {
    method2concordance = function(path, method) {
        score = scores[[method]]
        score[index$effect == "inhibiting",] = - score[index$effect == "inhibiting",]

        ranks = ar$map(score, along=2, function(x) {
            re = order(x)
            re[is.na(x)] = NA
            re
        })

        s1 = ranks[index$pathway == path, path]
        s2 = ranks[index$pathway != path, path]
        wilcox.test(s1, s2, paired=FALSE, alternative="greater") %>%
            broom::tidy() %>%
            as.list()
    }

    idx_df = df$create_index(path = unique(index$pathway),
                             method = names(scores),
                             expand_grid = TRUE)
    df$call(idx_df, function(path, method) method2concordance(path, method) %catch% NA) %>%
        mutate(adj.p = p.adjust(result, method="fdr"),
               logp = log(result))
}

hyper = ar$construct(result ~ path + method, data=test_hyper())
wilcox = test_wilcox()
effect = ar$construct(statistic ~ path + method, data=wilcox)
pval = ar$construct(p.value ~ path + method, data=wilcox)

pdf("validate.pdf", paper="a4r", width=26, height=20)
pheatmap::pheatmap(result,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE)
dev.off()
