library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
plt = import('plot')

# load scores
s1 = io$load('../../scores/gdsc/speed_norm.RData')
s2 = io$load('../../scores/gdsc/reactome.RData')
s3 = io$load('../../scores/gdsc/go.RData')

# load sanger data
Ys = gdsc$drug_response('IC50s') # or AUC
tissues = gdsc$tissues(minN=5)
ar$intersect(s1, s2, s3, tissues, Ys, along=1)

## save pdf w/ pan-cancer & tissue specific
pdf("VOLCANO.pdf", width=12, height=10)
on.exit(dev.off())

assocs2plot = function(scores, ylim) {
    assocs = st$lm(Ys ~ tissues + scores) %>%
        filter(term == "scores") %>%
        select(-term, -tissues) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               label = paste(Ys, scores, sep=":"))
    message(sum(assocs$adj.p < 0.05), " significant associations (5% FDR)")
    assocs %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
        plt$volcano(base.size=0.2, p=0.05, ceil=1e-10, xlim=c(-1,1), ylim=ylim)
}
a1 = assocs2plot(s1, c(1,1e-10))
a2 = assocs2plot(s2, c(1,8e-5))
a3 = assocs2plot(s3, c(1,8e-5))

left = arrangeGrob(a2, a3, ncol=1, nrow=2)
frame = arrangeGrob(left, a1, ncol=2, nrow=1)
print(frame)
