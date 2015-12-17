b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
plt = import('plot')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "plots/speed_matrix.pdf"
MUTFILE = "mutations_annotated_pathwayactivities_v3_mikeformat.txt"

scores = io$load(INFILE)
rownames(scores) = substr(rownames(scores), 1, 16)

mut = io$read_table(MUTFILE, header=TRUE) %>%
    transmute(hgnc = GENE_NAME,
              sample = substr(Tumor_Sample_Barcode, 1, 16),
              study = tcga$barcode2study(Tumor_Sample_Barcode)) %>%
    filter(!is.na(study) & study != "READ")

mut$mut = 1
mut = ar$construct(mut ~ sample + hgnc,
                   data=mut, fun.aggregate = length) > 0
ar$intersect(mut, scores)
mut = mut[,colSums(mut) > 35]
study = tcga$barcode2study(rownames(mut))

# associations
result = st$lm(scores ~ study + mut) %>%
    filter(term == "mutTRUE") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

# matrix plot
lims = result %>%
    filter(adj.p < 0.01) %>%
    select(estimate) %>% unlist() %>%
    abs() %>% max()

p1 = result %>%
    mutate(label = ifelse(p.value < 1e-3, "*", "")) %>%
    mutate(label = ifelse(p.value < 1e-20, "***", label)) %>%
    plt$cluster(estimate ~ scores + mut) %>%
    filter(adj.p < 0.1) %>%
    plt$matrix(estimate ~ scores + mut, color="estimate", limits=c(-lims,lims))

pdf(OUTFILE, paper="a4r", width=26, height=20)
print(p1)
dev.off()
