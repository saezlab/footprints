library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
cbp = import('data/cbioportal')

SCOREFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
VARFILE = commandArgs(TRUE)[2] %or% "gene_variants.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "assocs/speed_matrix.RData"

scores = io$load(SCOREFILE)
rownames(scores) = substr(rownames(scores), 1, 15) # do not include portion
scores = scores[!duplicated(rownames(scores)),]
gene_variants = io$load(VARFILE)

###
### do associations between scores and variants
###

# get samples where we have expression, mutation (and CNV?)
gene = as.data.frame(gene_variants$variants) # "duplicate row names" error otherwise!?
gene = gene[!duplicated(gene$id),] # this could be done better
ar$intersect(scores, gene$id)

# differences in pathway scores between studies could be anything (including
# the tissue just have the pathway more active irrespective of p53 status),
# so discard this and look only at effects where specific variants involved
variants = st$lm(scores ~ study * variant, data=gene) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    filter(! grepl("^study[^:]+$", term)) %>%
    arrange(adj.p) %>%
    mutate(label = paste(scores, sub("study", "", term), sep="_"))

###
### do associations between tissues and variants
###
scores = io$load(SCOREFILE)
rownames(scores) = substr(rownames(scores), 1, 15) # do not include portion
scores = scores[!duplicated(rownames(scores)),]

gene = ar$construct(variant ~ id + hgnc, data=gene_variants$cna, fill="normal",
    fun.aggregate=function(x) ifelse(length(unique(x) == 1), unique(x), "NA")) %>%
    as.data.frame()
for (cc in 1:ncol(gene))
    gene[[cc]] =  relevel(as.factor(gene[[cc]]), "normal")

study = setNames(gene_variants$cna$study, gene_variants$cna$id)
ar$intersect(scores, gene, study)

cnas = st$lm(scores ~ study * gene) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    filter(! grepl("^study[^:]+$", term)) %>%
    arrange(adj.p) %>%
    mutate(label = paste(scores, gene, sub("study", "", sub("gene", "", term)), sep="_"))

# save results to plot later
save(variants, cnas, file=OUTFILE)
