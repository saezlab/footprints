library(dplyr)
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')

edf =  io$read_table('speed2_norm_counts_hek293.tsv', header=TRUE) %>%
    select(-gene_id, -gene_type)
expr = data.matrix(edf[,2:ncol(edf)])
rownames(expr) = edf$gene_name

model = io$load('../../model/model_matrix.RData')$model
ar$intersect(model, expr, along=1)

scores = t(expr) %*% model %>%
    ar$map(along=1, scale)
rownames(scores) = sub("_[0-9]$", "", rownames(scores))

conditions = grep("BSA", unique(rownames(scores)), value=TRUE, invert=TRUE)
pathways = colnames(scores)

pval = matrix(NA, nrow=length(conditions), ncol=length(pathways),
              dimnames = list(conditions, pathways))
est = pval

for (cond in conditions)
    for (path in pathways) {
        ctl = paste0("BSA", sub("^[^_]+", "", cond))
        sctl = scores[rownames(scores) == ctl, path]
        sptb = scores[rownames(scores) == cond, path]
        tt = broom::tidy(t.test(sptb, sctl))
        pval[cond, path] = tt$p.value / 2 # one-sided test
        est[cond, path] = tt$estimate
    }

pval = ar$melt(pval) %>%
    transmute(ptb = Var1, pathway=Var2, p.value = value_df) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))
est = ar$melt(est) %>%
    transmute(ptb = Var1, pathway=Var2, estimate = value_df)
df = merge(pval, est, by=c('ptb','pathway'))

df %>%
    mutate(label = ifelse(p.value < 0.05, "*", "")) %>%
    plt$matrix(estimate ~ ptb + pathway, label=label)
