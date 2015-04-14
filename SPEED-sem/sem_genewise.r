#TODO: add sem code here
library(modules)
library(dplyr)
io = import('io')
ar = import('array')
ll = import('list')
hpc = import('hpc')
#an = import('anova')
#lm = import('lm')
#sg = import('sanger_robject')
#plt = import('plots')

# load zscore/dscore matrix for scores
# load index for pathways, convert to matrix
zscores = t(io$load('../SPEED-Data/SPEED2zmats.RData')$rma_none)
#dscores = t(io$load('../SPEED-Data/SPEED2dmats.RData')$rma_none)
index = io$read_table("../SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
    select(id, pathway, cells, GSE, GPL) %>%
    dplyr::filter(id %in% rownames(zscores))
zscores = zscores[index$id,]

# genes = matrix [samples x genes]
# pathways = character vector which pathway
createModel = function(genes, pathways) {
    fit1gene = function(gene, pathwayDf) {
        library(lavaan)
#        df = cbind(gene=gene[,1], pathwayDf) #TODO: ar$split should have drop=F/T
        df = cbind(gene=gene, pathwayDf) # works with hpc$Q (this is likely a "bug")
#FIXME: bsub error: gene incorrect dimensions
        model_str = paste(lapply(colnames(pathwayDf), function(p) {
            paste(p, "~", 'gene')
        }), collapse="\n")

#        print(model_str)
#        print(head(df))
        fit = suppressWarnings(sem(model_str, data=df)) # non-convergence->se=0 (?)

        idx = fit@ParTable$op=="~"
        est = fit@Fit@est[idx]
        se = fit@Fit@se[idx]
        name = fit@ParTable$lhs[idx]
        list(est = setNames(est, name),
             se = setNames(se, name))
    }

    # put both together in data.frame
    pathwayDf = as.data.frame(ar$mask(pathways))+0
    colnames(pathwayDf) = make.names(colnames(pathwayDf))

#    lapply(ar$split(zscores, along=2), function(g) fit1gene(g, pathwayDf))
    hpc$Q(fit1gene, gene=zscores, more.args=list(pathwayDf=pathwayDf),
          n.chunks=100, memory=1024)
# check: if I provide x-act as coeffs, should this work with lm as well?
    # construct the formulas according to scores+pathways
    #  for each grep ^path\\., set up regression:
    #  path ~ gene.1 + gene.2 + ...
    #  (this is standard multiple regression *only* to check results)
}

colnames(zscores) = make.names(colnames(zscores))
model = createModel(zscores, index$pathway)
result = ll$transpose(model)
save(result, file="sem_genewise_result.RData")

#TODO: create estimate/t-value table, select top genes, run assocs to compare to lm
