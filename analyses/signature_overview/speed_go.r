library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
gsea = import('../../util/gsea')
genesets = import('../../util/genesets')

#' Computes Gene Set Enrichment scores using piano
#'
#' @param col   Column index for 'mat'  
#' @param set   Name of the gene set
#' @param mat   Matrix of p-values, signed with the direction of change
#' @param sets  A list of character vectors for all gene sets
run_piano = function(col, set, mat, sets) {
    re = piano::runGSA(geneLevelStats = abs(mat[,col]),
                       direction = sign(mat[,col]),
                       gsc = piano::loadGSC(stack(sets[set])))

    sapply(re[grepl("^p[^A]", names(re))], c)
}

if (is.null(module_name())) {
    OUTFILE = commandArgs(TRUE)[1]

    assocs = io$load('../../model/model_matrix.RData')$assocs
    zfit = ar$construct(zscore ~ gene + pathway, data=assocs)
    pval = ar$construct(p.value ~ gene + pathway, data=assocs) * sign(zfit)

    go = genesets$get("go", mapped=FALSE, to_data_frame=FALSE)[["Gene Ontology"]]
    go = genesets$filter(go, valid=rownames(zfit), min=5, max=200)

    index = b$expand_grid(pathway = colnames(pval), set=names(go))
    re = clustermq::Q(run_piano, col=index$pathway, set=index$set,
                      const = list(mat=na.omit(pval), sets=go),
                      n_jobs = 100)
    index = cbind(index, ar$stack(re, along=1))

    save(index, file="speed_go.RData")
}
