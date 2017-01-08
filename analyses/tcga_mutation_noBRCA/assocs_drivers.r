# get pathway scores and mutations, and correlate them with each other
b = import('base')
io = import('io')
st = import('stats')
ar = import('array')
plt = import('plot')
tcga = import('data/tcga')
gdsc = import('data/gdsc')

subs2assocs = function(subs, mut, scores) {
    gdsc = import('data/gdsc')

    message(subs)
    if (grepl("pan", subs)) {
        drivers = unique(gdsc$drivers()$HGNC)
        m = filter(mut, Hugo_Symbol %in% drivers)
    } else {
        drivers = unique(gdsc$drivers(subs)$HGNC)
        m = filter(mut, Study==subs & Hugo_Symbol %in% drivers)
    }

    num_sample = length(unique(m$sample))
    altered = m$hgnc
    m$mut = 1
    m = ar$construct(mut ~ Tumor_Sample_Barcode + Hugo_Symbol,
                     data=m, fun.aggregate = length) > 0
    ar$intersect(m, scores)

    if (grepl("cov", subs)) {
        study = tcga$barcode2study(rownames(scores))
        assocs = st$lm(scores ~ study + m)
    } else {
        assocs = st$lm(scores ~ m)
    }

    result = assocs %>%
        filter(term == "mTRUE") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               label = paste(m, scores, sep=":"),
               subset = subs)
}

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/pathways_mapped/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "snp_drivers.pdf"

studies = setdiff(import('../../config')$tcga$tissues_with_normals, "BRCA")
driver_studies = unique(gdsc$drivers()$tissue)
scores = tcga$map_id(io$load(INFILE), along=1, id_type="specimen", subset="primary")

# AAChange is not avail in eg. BRCA (and others)
mut = tcga$mutations(id_type="specimen", subset="primary") %>%
    filter(Study %in% intersect(studies, driver_studies) &
           Tumor_Sample_Barcode %in% rownames(scores) &
           Variant_Classification != "Silent")

print(table(mut$Study))

methods = c("pan", "pan_cov", sort(unique(mut$Study)))
assocs = bind_rows(lapply(methods, function(x) subs2assocs(x, mut, scores)))

save(assocs, file=OUTFILE)
