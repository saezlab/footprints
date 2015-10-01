library(circlize)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/speed/gsea_reactome.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "perturb2expr.pdf"

# scale the same way as pheatmap does
scale_pheatmap = function(mat, n=100, center=TRUE) {
    color = 1:n
#    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(n)

    # scale matrix
    x = ar$map(mat, along=1, scale)

    # generate breaks
    m = max(abs(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))))
    breaks = seq(-m, m, length.out = n + 1)

    scaled = color[as.numeric(cut(c(x), breaks=breaks, include.lowest=TRUE))]
    mat[] = scaled
    mat
}

# load data
data = io$load(INFILE)
index = data$index
scores = scale_pheatmap(data$scores)
sign = ifelse(index$effect == "activating", 1, -1)
pathway = sign * t(ar$mask(index$pathway)) + 0

# ...
scores = scores/100 - 0.5

# MAPK-PI3K :: MAPK.E-GEOD.45757
##hist(scores[index$pathway == "MAPK",'PI3K'])
#blocked = scores[index$pathway == "MAPK",'PI3K'] < -0.2
#blocked = names(blocked)[blocked]
#as.data.frame(index[index$id %in% blocked,])
MAPK_PI3K = index[index$accession == "E-GEOD-45757",]
MAPK_PI3K$PI3K = scores[MAPK_PI3K$id,'PI3K']
MAPK_PI3K = as.data.frame(MAPK_PI3K)

# MAPK-Trail :: MAPK.E-GEOD-35230
##hist(scores[index$pathway == "MAPK",'Trail'])
#blocked = scores[index$pathway == "MAPK",'Trail'] < -0.05
#blocked = names(blocked)[blocked]
#as.data.frame(index[index$id %in% blocked,])
MAPK_Trail = index[index$accession == "E-GEOD-35230",]
MAPK_Trail$Trail = scores[MAPK_Trail$id,'PI3K']
MAPK_Trail = as.data.frame(MAPK_Trail)

# EGFR-TNFa

# VEGF-JAK-STAT

# compute associations
result = st$lm(scores ~ pathway) %>%
    filter(p.value < 0.05) %>%
    transmute(from=pathway, to=scores, value=effect)

chordDiagram(result, reduce=0.1)
