# point of this file:
# - use the zscores to create a linear model
library(dplyr)
io = import('io')
ar = import('array')
st = import('stats')
df = import('data_frame')
oc = import('./optimize_cutoff')

EXPR = commandArgs(TRUE)[1] %or% "../data/expr.RData"
ZSCORES = commandArgs(TRUE)[2] %or% "zscores.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "model_linear.RData"

data = io$load('../data/expr.RData')
records = data$records
expr = data$expr
zscores = io$load('zscores.RData')

index = records %>%
    lapply(function(x) x[sapply(x, length) == 1]) %>%
    bind_rows() %>%
    filter(is.na(exclusion))

# subset expr and records to train set
records = records[index$id]
expr = expr[index$id]

# prepare for model building
zscores = zscores[,index$id]
inh = index$effect=="inhibiting"
zscores[,inh] = -zscores[,inh]
zscores = t(zscores)

# fit model to pathway perturbations
mod = st$lm(zscores ~ 0 + pathway, data=index, min_pts=100,
            hpc_args=list(n_jobs=10, memory=2048)) %>%
    transmute(gene = zscores,
              pathway = sub("^pathway", "", term),
              zscore = estimate,
              p.value = p.value) %>%
    group_by(gene) %>%
    mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
    ungroup()

zfit = ar$construct(zscore ~ gene + pathway, data=mod)
pval = ar$construct(p.adj ~ gene + pathway, data=mod)

# take record, expr and calculate pathway scores
cut_df = df$create_index(zcut=seq(-2,2,0.1), pcut=c(1e-5, 1e-4, 1e-3, 1e-2, seq(0.05,0.5,0.05)),
                         args=list(zfit=zfit, pval=pval, records=records, expr=expr),
                         expand_grid=TRUE)
result = df$call(cut_df, oc$calc_concordance, hpc_args=list(n_jobs=50, memory=2048))

pdf("model_linear.pdf")
oc$plot_concordance(result)
dev.off()

# subset model matrix according to cutoffs
zfit[abs(zfit) < 1] = 0
zfit[pval > 0.15] = 0

# save resulting object
save(zfit, file=OUTFILE)
