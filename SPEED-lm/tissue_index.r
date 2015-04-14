library(modules)
b = import('base')
io = import('io')
ar = import('array')
sg = import('sanger_robject')
ll = import('list')
ent = import('./tissue_EN')

INFILE_MOD = commandArgs(TRUE)[1] %or% "mod_pathways-adrivers.RData"
INFILE_BG = commandArgs(TRUE)[2] %or% "mbg_pathways-adrivers.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "index_pathways-adrivers.RData"

ac = function(x) as.character(x)

# compose rsq extraction with load otherwise mem up
mod = io$load(INFILE_MOD)
rsq.real = ent$rsq2array(mod)
names(dimnames(rsq.real)) = c("drug","nfeat","tissue") #FIXME: do this in en_tissue
rsq.random = io$load(INFILE_BG)

calcPval = function(drug, nfeat, tissue) {
    vec2pval = function(vec, obs) {
        obs %or% return(1)
        fracNA = sum(is.na(vec))/length(vec)
        vec = log(b$omit$na(vec))
        pnorm(log(obs), mean=mean(vec), sd=sd(vec), lower.tail=F) * fracNA # empirical p-value
    }

    vec = rsq.random[drug, nfeat, tissue,]
    vec2pval(vec, rsq.real[drug, nfeat, tissue])
}

dn = dimnames(rsq.real)
pval = array(NA, dim=dim(rsq.real), dimnames=dn)

for (drug in dn[[1]])
    for (nfeat in dn[[2]])
        for (tissue in dn[[3]])
            pval[drug, nfeat, tissue] = calcPval(drug, nfeat, tissue)

pval = ar$map(pval, along=2, function(x) {
    re = rep(Inf, length(x))
    re[which(x==min(x))] = min(x)
    re
})
names(dimnames(pval)) = c("drug","nfeat","tissue") #FIXME: array should keep this!

# want: data.frame w/ tissue-drug-pval-nfeat + df w/ vectors
index = melt(pval)
index$drug = as.character(index$drug) #TODO: find a way to melt() without stringsAsFactors
index$tissue = as.character(index$tissue)

colnames(index) = c(names(dimnames(pval)), "pvalue")
index = arrange(dplyr::filter(index, pvalue<1), pvalue)

feats = sapply(1:nrow(index), function(i) {
    w = as.matrix(mod[[1]][[
        index$tissue[i]]][[index$drug[i]]]$weights)[,ac(index$nfeat[i])
    ] %catch% NA # when min CV error is not the same as min p-value
    w[w!=0]
})
rsquared = sapply(1:nrow(index), function(i) {
    rsq.real[index$drug[i], index$nfeat[i], index$tissue[i]] #TODO: find a better way to index this
})

index = cbind(index, rsquared=rsquared,
              feats=sapply(feats, f -> paste(names(f), collapse=",")))


# get additional drug data
Ys = sg$getDrugResponseForCellLines('IC50s') # or AUC
tissues = sg$getTissues(minN=5)
ar$intersect(tissues, Ys, along=1)

# add total variance for tissue-drug pair
# add min, max, median drug response for tissue-drug pair
add.feats = lapply(1:nrow(index), function(i) {
    resp = Ys[tissues==ac(index$tissue[i]), ac(index$drug[i])]
    list(
        resp_min = min(resp, na.rm=T),
        resp_max = max(resp, na.rm=T),
        resp_median = median(resp, na.rm=T),
        resp_var = var(resp, na.rm=T)
    )
})
add = ll$transpose(add.feats)

index$resp_min = add$resp_min
index$resp_max = add$resp_max
index$resp_median = add$resp_median
index$resp_var = add$resp_var

save(index, feats, file=OUTFILE)

