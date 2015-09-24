# point of this file:
# - use the zscores to create a linear model
library(dplyr)
b = import('base', attach_operators=FALSE)
import('base/operators')
io = import('io')
ar = import('array')
st = import('stats')

INFILE = commandArgs(TRUE)[1] %or% "../data/zscores.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "model_z0.RData"

# load speed data, index; filter for train set only
zobj = io$load('../data/scores.RData')
index = zobj$index %>% filter(is.na(exclusion))
zscores = zobj$zscores[,index$id]

# adjust object for linear modelling
inh = index$effect=="inhibiting"
zscores[,inh] = -zscores[,inh]
zscores = t(zscores)

# filter by mean abs z-score
meanz = ar$map(abs(zscores), along=2, function(x) {
    x = x[!is.na(x)]
    sum(x) / length(x)
})
library(mclust)
mc = Mclust(meanz)
keep = mc$classification == 1
index = index[keep,]
zscores = zscores[keep,]

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

# filter zfit to only include top 100 genes per pathway
zfit[apply(pval, 2, function(p) !b$min_mask(p, 100))] = 0

# save resulting object
save(zfit, file=OUTFILE)
