# get pathway scores and mutations, and correlate them with each other
b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
plt = import('plot')
tcga = import('data/tcga')

subs2plots = function(subs, mut, scores) {
    message(subs)
    if (grepl("pan", subs)) {
        m = mut %>%
            group_by(hgnc) %>%
            filter(n() >= 200) %>%
            ungroup()
        size = 0.5
    } else {
        m = filter(mut, study==subs) %>%
            group_by(hgnc) %>%
            filter(n() >= 10) %>%
            ungroup()
        size = 5
    }

    num_sample = length(unique(m$sample))
    altered = m$hgnc
    m$mut = 1
    m = ar$construct(mut ~ sample + hgnc,
                     data=m, fun.aggregate = length) > 0
    ar$intersect(m, scores)

    if (nrow(m) == 0) {
        warning("no overlap between mutations and scores for ", subs)
        return(NULL)
    }

    # associations
    if (grepl("cov", pan)) {
        study = tcga$barcode2study(rownames(scores))
        assocs = st$lm(scores ~ study + m)
    } else
        assocs = st$lm(scores ~ m)

    assocs %>%
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
        plt$volcano(base.size=size, p=0.1) +
            ggtitle(paste0(subs, " (", num_sample, " samples, ",
                           length(altered), " SNVs in ",
                           length(unique(altered)), " genes)"))
}

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "snp_all.pdf"

scores = io$load(INFILE)
rownames(scores) = substr(rownames(scores), 1, 16)

mut = tcga$mutations() %>%
    transmute(hgnc = Hugo_Symbol,
              sample = substr(Tumor_Sample_Barcode, 1, 16),
              study = study) %>%
    filter(!is.na(study) & study != "OV")

plots = mut$study %>%
    unique() %>%
    sort() %>%
    c("pan", "pan_cov", .) %>%
    lapply(function(s) subs2plots(s, mut, scores))

pdf(OUTFILE, paper="a4r", width=26, height=20)
for (plot in plots)
    print(plot)
dev.off()
