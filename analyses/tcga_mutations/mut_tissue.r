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
scores = io$load(INFILE) #TODO: add tissue normals with 0 mutations?
rownames(scores) = substr(rownames(scores), 1, 16)

mut = tcga$mutations() %>%
    transmute(hgnc = Hugo_Symbol,
              sample = substr(Tumor_Sample_Barcode, 1, 16),
              study = study) %>%
    group_by(study, hgnc) %>%
    filter(n() >= 5) %>%
    ungroup()

pdf("mut_tissue.pdf", paper="a4r", width=26, height=20)
on.exit(dev.off)

subs2plots = function(subs, mut, scores) {
    print(subs)
    m = filter(mut, study==subs)
    m$mut = 1
    m = ar$construct(mut ~ sample + hgnc, data=m, fun.aggregate = length) > 0
    ar$intersect(m, scores)

    # associations
    result = st$lm(scores ~ m) %>%
        filter(term == "mTRUE") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))

    # matrix plot
    p1 = result %>%
        mutate(label = ifelse(adj.p < 0.01, "*", "")) %>%
        plt$cluster(estimate ~ scores + m) %>%
        filter(adj.p < 0.1) %>%
        plt$matrix(estimate ~ scores + m, color="estimate") + ggtitle(subs)
    print(p1)

    # volcano plot
    p2 = result %>%
        mutate(label = paste(m, scores, sep=":")) %>%
        plt$color$p_effect(pvalue = "adj.p") %>%
        plt$volcano(base.size=0.1) + ggtitle(subs)
    print(p2)
}
lapply(unique(mut$study), function(s) subs2plots(s, mut, scores) %catch% NULL)
