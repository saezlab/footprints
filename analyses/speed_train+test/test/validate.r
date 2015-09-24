library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
df = import('data_frame')
plt = import('plot')

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
        
        list(p.value = 1.0 - phyper(q = sum(concordance) - 1,
                                    m = length(concordance),
                                    n = length(concordance) * ncol(score),
                                    k = length(concordance)),
             fraction = sum(concordance) / length(concordance))
    }

    idx_df = df$create_index(path = unique(index$pathway),
                             method = names(scores),
                             expand_grid = TRUE)
    df$call(idx_df, method2concordance) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))
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

hyper = test_hyper() %>%
    mutate(fraction = ifelse(p.value < 0.2, fraction, NA)) %>%
    mutate(label = ifelse(p.value < 0.05, "*", "")) %>%
    mutate(label = ifelse(p.value < 1e-5, "***", label)) %>%
    plt$matrix(fraction ~ path + method, palette="Blues")

pdf("validate.pdf", paper="a4r", width=26, height=20)
print(hyper)
dev.off()
