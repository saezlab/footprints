b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
plt = import('plot')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/pathways_mapped/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "heatmap_speed_matrix.pdf"

scores = io$load(INFILE)
rownames(scores) = substr(rownames(scores), 1, 16)
index = tcga$barcode2index(rownames(scores)) %>%
    filter(Short.Letter.Code %in% c("NT", "TP")) %>%
    mutate(tumor = Short.Letter.Code == "TP") %>%
    filter(!Study.Abbreviation %in% c("KICH", "KIRP")) %>%
    group_by(Study.Abbreviation) %>%
    filter(sum(!tumor) >= 10) %>%
    ungroup()

scores = scores[index$Bio.ID,]

# associations
result = st$lm(scores ~ tumor, data=index, subsets=index$Study.Abbreviation) %>%
    filter(term == "tumorTRUE") %>%
    select(-term) %>%
    na.omit() %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

# matrix plot
lims = result %>%
    filter(adj.p < 0.01) %>%
    select(estimate) %>% unlist() %>%
    abs() %>% max()

p1 = result %>%
    mutate(label = ifelse(p.value < 1e-3, "*", "")) %>%
    mutate(label = ifelse(p.value < 1e-20, "***", label)) %>%
    plt$cluster(estimate ~ scores + subset) %>%
    filter(adj.p < 0.1) %>%
    plt$matrix(estimate ~ scores + subset, color="estimate",
               limits=c(-lims,lims), palette="RdBu", color_rev=TRUE)

pdf(OUTFILE, paper="a4r", width=26, height=20)
print(p1)

library(pheatmap)
mat = ar$construct(estimate ~ scores + subset, data=result)
pheatmap(mat)

dev.off()
