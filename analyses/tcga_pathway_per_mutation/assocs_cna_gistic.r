# get pathway scores and mutations, and correlate them with each other
b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
plt = import('plot')
tcga = import('data/tcga')
gdsc = import('data/gdsc')

subs2assocs = function(subs, cna, scores) {
    message(subs)
    if (grepl("pan", subs)) {
        m = cna %>%
            group_by(hgnc) %>%
            filter(n() >= 50) %>%
            ungroup()
        size = 0.1
    } else {
        m = filter(cna, study==subs) %>%
            group_by(hgnc) %>%
            filter(n() >= 5) %>%
            ungroup()
        size = 0.5
    }

    num_sample = length(unique(m$barcode))
    m = ar$construct(gistic ~ barcode + hgnc, data=m,
                     fun.aggregate = mean, fill=0)
    ar$intersect(m, scores)

    if (nrow(m) == 0) {
        warning("no overlap between mutations and scores for ", subs)
        return(NULL)
    }

    # associations
    if (grepl("cov", subs)) {
        study = tcga$barcode2study(rownames(scores))
        assocs = st$lm(scores ~ study + m)
    } else
        assocs = st$lm(scores ~ m)

    result = assocs %>%
        filter(term == "m") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               label = paste(m, scores, sep=":"),
               subset = subs)
}

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/pathways_mapped/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "cna_gistic.RData"

scores = io$load(INFILE)
rownames(scores) = substr(rownames(scores), 1, 15)
cna = tcga$cna() %>%
    select(barcode, hgnc, gistic) %>%
    mutate(barcode = substr(barcode, 1, 15)) %>%
    filter(barcode %in% rownames(scores) & hgnc %in% gdsc$drivers()$HGNC) %>%
    unique() %>%
    mutate(study = tcga$barcode2study(barcode))

methods = c("pan", "pan_cov", sort(unique(cna$study)))
assocs = bind_rows(lapply(methods, function(x) subs2assocs(x, cna, scores)))

save(assocs, file=OUTFILE)
