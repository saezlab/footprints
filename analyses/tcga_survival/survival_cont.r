library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_norm.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "cont_speed_norm.pdf"

clinical = tcga$clinical() %>%
    transmute(study = study,
              age_days = - as.integer(patient.days_to_birth),
              alive = 1 - as.integer(is.na(patient.days_to_death)),
              surv_days = as.integer(patient.days_to_death %or%
                                     patient.days_to_last_followup),
              barcode = toupper(patient.bcr_patient_barcode),
              gender = as.factor(patient.gender)) %>%
    filter(surv_days > 0) #FIXME:

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

scores = io$load(INFILE)
# subset primary tumors only
#TODO: what if i do average over all tumor samples for each patient
#  (but strong signals should be in there either way)
scores = scores[substr(rownames(scores), 14, 16) == "01A",]
# construct clinical df for samples
cc = do.call(rbind, lapply(rownames(scores), function(s) {
    row = clinical[substr(s, 1, 12) == clinical$barcode,]
    if (nrow(row) == 1)
        row
    else
        NA
}))
nas = is.na(cc$study)
clinical = cc[!nas,]
scores = scores[!nas,]

if (nrow(clinical) < 10)
    stop("survival+expression scores < 10 observations, stopping")

clinical = as.list(clinical)
clinical$scores = scores

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

# tissue covariate
#TODO: add gender as covar; but: util tries to subset it, shouldn't
st$coxph(surv_days + alive ~ age_days + study + scores, data=clinical, min_pts=100) %>%
    filter(term == "scores") %>%
    select(scores, estimate, p.value, size) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    plt$color$p_effect("adj.p") %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1) %>%
    print()

# separate regressions for each tissue
#TODO: add gender + make it work w/ only one
st$coxph(surv_days + alive ~ age_days + scores, subsets=clinical$study, data=clinical, min_pts=20) %>%
    filter(term == "scores") %>%
    select(scores, subset, estimate, p.value, size) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup() %>%
    plt$color$p_effect("adj.p") %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) + ggtitle(sum(clinical$adj.p < 0.1)) %>%
    print()
