# point of this file:
# - use the zscores to create a linear model
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')

INFILE = commandArgs(TRUE)[1] %or% "../data/zscores.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "model_iter.RData"

# load speed data, index
zobj = io$load('../data/zscores.RData')
zscores = zobj$scores
index = zobj$index
inh = index$effect=="inhibiting"
zscores[,inh] = -zscores[,inh]
zscores = t(zscores)

# load dscores
dobj = io$load('../data/dscores.RData')
dscores = dobj$scores
dscores[,inh] = -dscores[,inh]
dscores = t(dscores)

iterate = function(zscores, dscores, index) {
    # fit model to pathway perturbations
    data = c(as.list(index), list(zscores=zscores))
    mod = st$lm(zscores ~ 0 + pathway, data=data) %>%
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
    zfit[apply(pval, 2, p -> !b$min_mask(p, 100))] = 0

    # make sure scores are highest in perturbed pathway, discard otherwise
    scores = dscores[,rownames(zfit)] %*% zfit
    bscores = ar$map(scores, along=2, x -> x==max(x))
    bscores[,"EGFR"] = bscores[,"MAPK"] = bscores[,"EGFR"] | bscores[,"MAPK"]
    keep = as.logical(rowSums(bscores & ar$mask(index$pathway)))
    print(sum(keep))

    if (all(keep))
        zfit
    else
        iterate(zscores[keep,], dscores[keep,], index[keep,])
}

# build linear model iteratively discarding bad arrays
zfit = iterate(zscores, dscores, index)

# save resulting object
save(zfit, file=OUTFILE)
