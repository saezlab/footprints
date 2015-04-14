library(modules)
io = import('io')
gs = import('genesets')
gsea = import('gsea')
icgc = import('icgc')

sets = gs$list(
    GO = gs$list('Ontologies')$GO_Biological_Process
    BioCarta = gs$list('Pathways')$BioCarta
    gatza = gs$list('gatza')$gatza
    TFs = gs$list('TF_iorio')
)

# load mutational data
mut = icgc$getMutations()
# (1) known mutations (GO categories)
# (2) pathway mutations

# load expression data
expr = icgc$getRNASeq()
# (3) gatza GSEA
# (4) TF regulons
# (5) GO categories (expression)
# (6) speed lm
# (7) speed sem
# load RPPA data
# (8) match them up to the RPPA profiles (calc correlation coeffs)

#NEED: TCGA mutational data; pathway<->gene mapping
