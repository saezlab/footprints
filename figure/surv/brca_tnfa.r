library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
icgc = import('data/icgc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "BRCA.pdf"

#    mutate(known_survival_time = donor_survival_time %or% donor_interval_of_last_followup) %>%
clinical = icgc$clinical() %>%
    filter(#!is.na(donor_survival_time),
           donor_survival_time != 0,
           donor_vital_status %in% c('alive','deceased'),
           grepl("Primary tumour", specimen_type),
           grepl("Breast.*US", study)) %>%
    select(icgc_specimen_id, donor_survival_time, donor_vital_status,
           donor_sex, donor_age_at_diagnosis)

scores = io$load(INFILE)
rownames(clinical) = clinical$icgc_specimen_id
ar$intersect(scores, clinical, along=1)

scores = ar$map(scores, along=1, subsets=clinical$tissue, function(x) {
    x = st$median_scale(x)
    re = ifelse(x > 1, TRUE, FALSE)
    if (sum(re) < 5)
        rep(NA, length(x))
    else
        re
})

if (nrow(clinical) < 10)
    stop("survival+expression scores < 10 observations, stopping")

clinical = as.list(clinical)
clinical$status = as.integer(clinical$donor_vital_status == "alive")
clinical$time = clinical$known_survival_time
clinical$scores = scores

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

st$coxph(time + status ~ scores, data=clinical, min_pts=20) %>%
    select(-time, -status, -term) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup() %>%
    plt$color$p_effect() %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) %>%
    print()
