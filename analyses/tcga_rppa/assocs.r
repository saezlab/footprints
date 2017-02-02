b = import('base')
io = import('io')
st = import('stats')
tcga = import('data/tcga')

INFILE = "../../scores/tcga/pathways_mapped/speed_matrix.RData"
OUTFILE = "assocs.RData"

# load expression, RPPA for all cancers where both available
rppa = t(tcga$rppa())
scores = io$load(INFILE)
tcga$intersect(rppa, scores, along=1)

tissues = tcga$barcode2study(rownames(scores))

result = st$lm(scores ~ rppa, subsets=tissues, hpc_args=list(n_jobs=10)) %>%
    filter(term == "rppa") %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

save(result, file=OUTFILE)
