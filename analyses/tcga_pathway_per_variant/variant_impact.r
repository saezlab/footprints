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
gene_variants = io$load(VARFILE)

###
### do associations between scores and variants
###

# get samples where we have expression, mutation (and CNV?)
gene = gene_variants$variants %>%
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
variants = st$lm(scores ~ study * variant, data=gene) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    filter(! grepl("^study[^:]+$", term)) %>%
    arrange(adj.p) %>%
    mutate(label = paste(scores, sub("study", "", term), sep="_")) %>%
    plt$color$p_effect(pvalue="adj.p", thresh=0.1)

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
    mutate(label = paste(scores, gene, sub("study", "", sub("gene", "", term)), sep="_")) %>%
    plt$color$p_effect(pvalue="adj.p", thresh=0.1)

###
### volcano plots of the results
###

pdf(OUTFILE)#, width=26, height=20)

variants %>%
    filter(!is.na(size)) %>%
    plt$volcano(base.size=10, p=0.1) + ggtitle("pan-cancer variants")

variants %>%
    filter(is.na(size)) %>%
    mutate(size = 50) %>%
    plt$volcano(p=0.1) + ggtitle("tissue-specific variants (interaction terms)")

cnas %>%
    filter(!is.na(size)) %>%
    plt$volcano(base.size=0.01, p=0.1) + ggtitle("pan-cancer CNAs")

cnas %>%
    filter(is.na(size)) %>%
    mutate(size = 50) %>%
    plt$volcano(p=0.1) + ggtitle("tissue-specific CNAs (interaction terms)")

dev.off()
