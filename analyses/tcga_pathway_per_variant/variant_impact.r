library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
cbp = import('data/cbioportal')
plt = import('plot')

SCOREFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
VARFILE = commandArgs(TRUE)[2] %or% "gene_variants.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "speed_matrix.pdf"

scores = io$load(SCOREFILE)
rownames(scores) = substr(rownames(scores), 1, 15) # do not include portion
scores = scores[!duplicated(rownames(scores)),]
#                c('p53','PI3K','MAPK','Hypoxia','NFkB')]

# get samples where we have expression, mutation (and CNV?)
gene = io$load(VARFILE)$variants %>%
#    filter(study != "BRCA") %>% # to set if TP53_R282W w/o BRCA: p<0.0012 (4% FDR)
#    filter(hgnc == "TP53") %>%
    mutate(variant = paste(hgnc, variant, sep="_")) %>%
    group_by(variant) %>%
    mutate(n=n()) %>%
    ungroup() %>%
    mutate(variant = ifelse(n>=20 & !grepl("MUTATED", variant),
                            variant, "other")) %>%
    mutate(variant = relevel(as.factor(variant), "other")) %>%
    as.data.frame() # error otherwise; why?

gene = gene[!duplicated(gene$id),] # this could be done better
ar$intersect(scores, gene$id)

# differences in pathway scores between studies could be anything (including
# the tissue just have the pathway more active irrespective of p53 status),
# so discard this and look only at effects where specific variants involved
result = st$lm(scores ~ study * variant, data=gene) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    filter(! grepl("^study[^:]+$", term)) %>%
    arrange(adj.p) %>%
    mutate(label = paste(scores, sub("study", "", term), sep="_")) %>%
    plt$color$p_effect(pvalue="adj.p", thresh=0.1)

pdf(OUTFILE)#, width=26, height=20)

result %>%
    filter(!is.na(size)) %>%
    plt$volcano(base.size=10, p=0.1) + ggtitle("variants")

result %>%
    filter(is.na(size)) %>%
    mutate(size = 50) %>%
    plt$volcano(p=0.1) + ggtitle("interactions")

dev.off()
