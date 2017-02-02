b = import('base')
io = import('io')
st = import('stats')
tcga = import('data/tcga')

INFILE = "../../scores/tcga/pathways_mapped/speed_matrix.RData"
OUTFILE = "assocs.RData"

# load required pathway scores and rppa data
rppa = t(tcga$rppa())
scores = io$load(INFILE)
tcga$intersect(rppa, scores, along=1)

tissues = tcga$barcode2study(rownames(scores))

# calculate associations across and within tissues
pan = st$lm(scores ~ tissues + rppa) %>%
    filter(term == "rppa") %>%
    select(-term)

tissue = st$lm(scores ~ rppa, subsets=tissues, hpc_args=list(n_jobs=10)) %>%
    filter(term == "rppa") %>%
    select(-term) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

save(pan, tissue, file=OUTFILE)
