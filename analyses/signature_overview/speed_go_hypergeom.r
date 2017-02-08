library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
gsea = import('../../util/gsea')
genesets = import('../../util/genesets')

OUTFILE = commandArgs(TRUE)[1] %or% "speed_go_hypergeom.RData"

# we're loading this twice here just to get ground truth
# could be more efficient ..
assocs = io$load('../../model/model_matrix.RData')$assocs
zfit = ar$construct(zscore ~ gene + pathway, data=assocs)

speed = genesets$get("go", mapped=TRUE, to_data_frame=FALSE)$PROGENy
go = genesets$get("go", mapped=FALSE, to_data_frame=FALSE)[["Gene Ontology"]]
go = genesets$filter(go, valid=rownames(zfit), min=5, max=200)

result = st$hypergeometric_test(rownames(zfit), go, speed)
result[result == 0] = .Machine$double.eps
save(result, file=OUTFILE)
