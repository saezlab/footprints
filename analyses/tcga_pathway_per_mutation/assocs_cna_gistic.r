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
    study = tcga$barcode2study(rownames(scores))

    if (grepl("pan", subs)) {
        m = ar$filter(cna, along=2, function(x) sum(x!=0) > 50, subsets=study)
        size = 0.1
    } else {
        m = ar$filter(cna, along=2, function(x) sum(x!=0) > 5, subsets=study)
        size = 0.5
    }

    ar$intersect(m, scores)

    if (nrow(m) == 0) {
        warning("no overlap between CNA and scores for ", subs)
        return(NULL)
    }

    # associations
    if (grepl("cov", subs)) {
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
CNAFILE = commandArgs(TRUE)[2] %or% "cna_driver_matrix.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "cna_gistic.RData"

scores = io$load(INFILE)
cna = io$load(CNAFILE)
studies = unique(tcga$barcode2study(rownames(cna)))

methods = c("pan", "pan_cov", sort(studies))
assocs = bind_rows(lapply(methods, function(x) subs2assocs(x, cna, scores)))

save(assocs, file=OUTFILE)
