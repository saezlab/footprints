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
    study = tcga$barcode2study(rownames(cna))

    if (grepl("pan", subs)) {
        m = ar$filter(cna, along=1, function(x) sum(x!=0, na.rm=TRUE) > 50, subsets=study)
        size = 0.1
    } else {
        m = ar$filter(cna[study==subs,], along=1, function(x) sum(x!=0, na.rm=TRUE) > 5)
        size = 0.5
    }

    if (nrow(m) == 0) {
        warning("no overlap between CNA and scores for ", subs)
        return(NULL)
    }

    ar$intersect(scores, m, along=1)
    study = tcga$barcode2study(rownames(m))

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

studies = setdiff(import('../../config')$tcga$tissues_with_normals, "BRCA")
scores = tcga$map_id(io$load(INFILE), along=1, id_type="specimen", subset="primary")
cna = io$load(CNAFILE)

methods = c("pan", "pan_cov", sort(studies))
assocs = bind_rows(lapply(methods, function(x) subs2assocs(x, cna, scores)))

save(assocs, file=OUTFILE)
