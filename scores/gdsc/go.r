# this needs
# * a gsea script
# * gene ontology categories + genes in them
# * perform gsea
# * save the scores for each cell line in robject
b = import('base')
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')

OUTFILE = commandArgs(TRUE)[2] %or% "go.RData"

# load sanger data
Ys = gdsc$getDrugResponse('IC50s') # or AUC
tissues = gdsc$getTissues(minN=5)
expr = t(gdsc$getBasalExpression())
ar$intersect(expr, tissues, Ys, along=1)

# calculate scores on panel
scores = expr %*% zfit #%>%
#    ar$map(along=1, base::scale) # pathway across experiments

save(scores, file=OUTFILE)
