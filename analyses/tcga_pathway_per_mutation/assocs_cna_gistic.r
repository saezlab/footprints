# get pathway scores and mutations, and correlate them with each other
b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
plt = import('plot')
tcga = import('data/tcga')

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

    num_sample = length(unique(m$sample))
    altered = m$hgnc
    m$altered = 1
    m = ar$construct(altered ~ sample + hgnc, data=m,
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

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "cna_gistic.RData"
CNAFILE = "cna.txt"

scores = io$load(INFILE)
rownames(scores) = substr(rownames(scores), 1, 15)
gistic_lookup = setNames(c("+", "-"), c(2, -2))

cna = io$read_table(CNAFILE, header=TRUE) %>%
    transmute(hgnc = GENE_NAME,
              sample = substr(Tumor_Sample_Barcode, 1, 15), # NO PORTION
              study = study,
              gistic = sapply(CNA_gistic, function(x) gistic_lookup[as.character(x)])) %>%
    mutate(hgnc = paste0(hgnc, gistic)) %>%
    select(-gistic) %>%
    filter(!study %in% c("KICH","LAML")) # no alteration present n>=cutoff

methods = c("pan", "pan_cov", sort(unique(cna$study)))
assocs = bind_rows(lapply(methods, function(x) subs2assocs(x, cna, scores)))

save(assocs, file=OUTFILE)
