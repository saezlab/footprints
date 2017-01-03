b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
plt = import('plot')
tcga = import('data/tcga')

INFILE = "../../scores/tcga/speed_linear.RData"
OUTFILE = "speed_linear.pdf"

# load expression, RPPA for all cancers where both available
cc = tcga$clinical() %>%
    transmute(code = toupper(patient.bcr_patient_barcode),
              study = study)
rppa = t(tcga$rppa())
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

pdf("protein_fits.pdf", paper="a4r", width=26, height=20)
on.exit(dev.off)

for (subs in unique(result$subset)) {
    print(subs)
    r2 = result %>% #TODO: cluster by all tissues, so order is the same?
        filter(subset == subs) %>%
        plt$cluster(estimate ~ scores + rppa) %>%
        group_by(rppa) %>%
        mutate(p.sum = sum(-log(adj.p))) %>%
        ungroup() %>%
        filter(adj.p < 0.05)

    keep = unique(r2[c('rppa','p.sum')]) %>%
        top_n(40, p.sum)

    r3 = r2 %>%
        filter(rppa %in% keep$rppa) %>%
        plt$matrix(estimate ~ scores + rppa, color="estimate") + ggtitle(subs)

    print(r3) %catch% NULL
}
