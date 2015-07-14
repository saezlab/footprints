library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
icgc = import('data/icgc')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "../../expr_cluster/speed_cluster.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "surv_speed.pdf"

# load scores
clusters = io$load(INFILE)
clusters = clusters[lapply(clusters, length) != 0]
scores = lapply(clusters, function(cl) {
    re = cl %>%
        lapply(function(x) x[grepl("^SP[0-9]+$", rownames(x)),]) %>%
        ar$stack(along=2)
    re[,colSums(re)>0]
})

# load sanger data
# separate associations for each tissue
clinical = icgc$clinical() %>% #TODO: select tumor only
#    filter(tumour_confirmed == "yes") %>%
#    filter(specimen_donor_treatment_type == "no treatment") %>%
    mutate(known_survival_time = donor_survival_time %or% donor_interval_of_last_followup) %>%
    select(icgc_specimen_id, known_survival_time, donor_vital_status,
           donor_sex, donor_age_at_diagnosis, tissue) %>%
    filter(!is.na(known_survival_time),
           donor_vital_status %in% c('alive','deceased'))

# save pdf w/ pan-cancer & tissue specific
pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

for (cur_tissue in names(scores)) {
    score = scores[[cur_tissue]]
    clin = filter(clinical, tissue == cur_tissue)
    rownames(clin) = clin$icgc_specimen_id
    ar$intersect(score, clin)

    if (nrow(clin) == 0)
        next

    assocs.tissue = clin %>%
        mutate(status = donor_vital_status=="alive",
               time = known_survival_time) %>%
        st$coxph(time + status ~ score, data=.) %>%
        select(-time, -status, -term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               label = score) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
        plt$volcano(p=0.1) + ggtitle(cur_tissue)
    print(assocs.tissue)
}
