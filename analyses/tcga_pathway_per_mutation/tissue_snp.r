# get pathway scores and mutations, and correlate them with each other
b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
plt = import('plot')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "tissue.pdf"
MUTFILE = "mutations_annotated_pathwayactivities_v3_mikeformat.txt"

scores = io$load(INFILE)
rownames(scores) = substr(rownames(scores), 1, 16)

mut = io$read_table(MUTFILE, header=TRUE) %>%
    transmute(hgnc = GENE_NAME,
              sample = substr(Tumor_Sample_Barcode, 1, 16),
              study = tcga$barcode2study(Tumor_Sample_Barcode)) %>%
    filter(!is.na(study)) %>%
    group_by(study, hgnc) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    filter(!study %in% c("READ"))

subs2plots = function(subs, mut, scores) {
    message(subs)
    m = filter(mut, study==subs)
    m$mut = 1
    m = ar$construct(mut ~ sample + hgnc,
                     data=m, fun.aggregate = length) > 0
    ar$intersect(m, scores)

    # associations
    result = st$lm(scores ~ m) %>%
        filter(term == "mTRUE") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))

#    # matrix plot
#    p1 = result %>%
#        mutate(label = ifelse(adj.p < 0.01, "*", "")) %>%
#        plt$cluster(estimate ~ scores + m) %>%
#        filter(adj.p < 0.1) %>%
#        plt$matrix(estimate ~ scores + m, color="estimate") +
#            ggtitle(subs)
#    print(p1)

    # volcano plot
    result %>%
        mutate(label = paste(m, scores, sep=":")) %>%
        plt$color$p_effect(pvalue="adj.p", thresh=0.1) %>%
        plt$volcano(base.size=0.1, p=0.1) + ggtitle(subs)
}

plots = mut$study %>%
    unique() %>%
    sort() %>%
    lapply(function(s) subs2plots(s, mut, scores))

pdf(OUTFILE, paper="a4r", width=26, height=20)
for (plot in plots)
    print(plot)
dev.off()
