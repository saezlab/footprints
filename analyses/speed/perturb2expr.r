library(circlize)
library(RColorBrewer)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/speed/gsea_reactome.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "perturb2expr.pdf"

file2assocs = function(fname) {
    # load data
    data = io$load(fname)
    index = data$index
    scores = data$scores
    sign = ifelse(index$effect == "activating", 1, -1)
    pathway = sign * t(ar$mask(index$pathway)) + 0

    # compute associations
    st$lm(scores ~ pathway)
}

reactome = file2assocs("../../scores/speed/gsea_reactome.RData") %>%
    filter(p.value < 1e-4) %>%
    transmute(from=pathway, to=scores, value=estimate)
speed = file2assocs("../../scores/speed/speed_matrix.RData") %>%
    filter(p.value < 1e-4) %>%
    transmute(from=pathway, to=scores, value=estimate)

all_pathways = unique(c(speed$from, speed$to, reactome$from, reactome$to))
color = setNames(brewer.pal(12, "Paired"), all_pathways)

pdf(OUTFILE, width=12, height=12)
chordDiagram(speed, grid.col=color, grid.border="#dedede", self.link=1, link.border="#dedede")
chordDiagram(reactome, grid.col=color, grid.border="#dedede", self.link=1, link.border="#dedede")
dev.off()
