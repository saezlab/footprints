library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
st = import('stats')
gdsc = import('data/gdsc')
plt = import('plot')

scores2assocs = function(scores) {
    assocs = st$lm(Ys ~ tissues + scores) %>%
        filter(term == "scores") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               label = paste(Ys, scores, sep=":"))
    message(sum(assocs$adj.p < 0.05), " significant associations (5% FDR)")
    assocs
}

conditional_assocs = function(scores) {
    best_speed_hit = a_speed %>%
        group_by(Ys) %>%
        filter(adj.p == min(adj.p)) %>%
        select(Ys, scores, adj.p) %>%
        distinct() # duplicate drug names in Ys

    condition_on = b$match(colnames(Ys), from=best_speed_hit$Ys, best_speed_hit$scores)
    condition = s_speed[,condition_on]
    colnames(condition) = names(condition_on)

    # return df with associations still significant after conditioning
    st$lm(Ys ~ tissues + condition + scores, group=c("Ys", "condition")) %>%
        filter(term == "scores") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        filter(adj.p < 0.05)
}

assocs2plot = function(assocs, ylim, condition=NULL) {
    assocs = assocs %>%
         plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1)

    if (!is.null(condition)) {
        subset_valid = condition %>%
            select(Ys, scores) %>%
            df$contains(assocs, .) %catch% FALSE

        assocs[!subset_valid & assocs$adj.p < 0.05,'color'] = "#00000000"
        assocs[!subset_valid & assocs$adj.p < 0.05,'circle'] = "#767676ff"
    }

    assocs %>%
        plt$volcano(base.size=0.2, p=0.05, ceil=1e-10, xlim=c(-1,1), ylim=ylim)
}

# load scores; s_: prefix for scores
s_speed = io$load('../../scores/gdsc/speed_matrix.RData')
s_reac = io$load('../../scores/gdsc/reactome.RData')
s_go = io$load('../../scores/gdsc/go.RData')

# load sanger data
Ys = gdsc$drug_response('IC50s') # or AUC
tissues = gdsc$tissues(minN=5)
ar$intersect(s_speed, s_reac, s_go, tissues, Ys, along=1)

# a_: prefix for associations
a_speed = scores2assocs(s_speed) # limit are there so axes are equal
a_reac = scores2assocs(s_reac)
a_go = scores2assocs(s_go)

# c_: prefix conditional associations
c_reac = conditional_assocs(s_reac)
c_go = conditional_assocs(s_go)

# v_: prefix for volcano plots
v_speed = assocs2plot(a_speed, c(0.57,1e-10)) # limit are there so axes are equal
v_reac = assocs2plot(a_reac, c(1,3e-5), condition=c_reac)
v_go = assocs2plot(a_go, c(1,3e-5), condition=c_go)

# arrange plots
right = arrangeGrob(v_go, v_reac, ncol=1, nrow=2)
frame = arrangeGrob(v_speed, right, ncol=2, nrow=1)

# save to pdf
pdf("conditional.pdf", width=12, height=10)
plot(frame)
dev.off()
