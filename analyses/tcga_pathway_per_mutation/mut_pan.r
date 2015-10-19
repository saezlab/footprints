# do the same as with the rppa data:
# get pathway scores and mutations, and correlate them with each other
b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
plt = import('plot')
tcga = import('data/tcga')

INFILE = "../../scores/tcga/speed_linear.RData"
OUTFILE = "speed_linear.pdf"

# load expression, RPPA for all cancers where both available
mut = tcga$mutations() %>%
    transmute(hgnc = Hugo_Symbol,
              sample = substr(Tumor_Sample_Barcode, 1, 16),
              study = study) %>%
    group_by(study, hgnc) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    ar$construct(1 ~ hgnc + sample, data=mut)

scores = io$load(INFILE)
rownames(rppa) = substr(rownames(rppa), 1, 16)
rownames(scores) = substr(rownames(scores), 1, 16)
ar$intersect(rppa, scores)
tissues = cc$study[match(substr(rownames(rppa), 1, 12), cc$code)]

result = st$lm(scores ~ rppa, subsets=tissues) %>%
    filter(term == "rppa") %>%
    select(-term, -std.error, -statistic) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()
