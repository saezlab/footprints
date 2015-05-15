library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
icgc = import('data/icgc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

# load clinical data, extract survival time+status
clinicalsample = icgc$clinical_sample() %>%
    select(icgc_specimen_id, icgc_sample_id)

clinical = icgc$clinical() %>% #TODO: select tumor only
#    filter(tumour_confirmed == "yes") %>%
#    filter(specimen_donor_treatment_type == "no treatment") %>%
    mutate(known_survival_time = donor_survival_time %or% donor_interval_of_last_followup) %>%
    select(icgc_specimen_id, known_survival_time, donor_vital_status,
           donor_sex, donor_age_at_diagnosis, tissue) %>%
    filter(!is.na(known_survival_time),
           donor_vital_status %in% c('alive','deceased')) %>%
    left_join(clinicalsample, by="icgc_specimen_id") %>%
    filter(icgc_sample_id %in% icgc$names$rna_seq()[[1]])

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - filter for last known alive? [where is this field?]
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

scores = io$load(INFILE)
survival = clinical[c('known_survival_time','donor_vital_status','tissue')]
rownames(survival) = clinical$icgc_sample_id
survival$donor_vital_status = as.integer(survival$donor_vital_status=="alive")
survival = as.matrix(survival)
ar$intersect(scores, survival, along=1)
status = as.integer(survival[,'donor_vital_status'])
time = as.integer(survival[,'known_survival_time'])
tissue = survival[,'tissue']

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

# tissue covariate
st$coxph(time + status ~ tissue + scores) %>%
    filter(term == "scores") %>%
    select(-time, -status, -tissue, -term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    plt$color$p_effect() %>%
    mutate(label = scores) %>%
    plt$volcano() %>%
    print()

# run cox regression for associations (tissue as covariate)
st$coxph(time + status ~ scores, subsets=tissue) %>%
    select(-time, -status, -term) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup() %>%
    plt$color$p_effect() %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano() %>%
    print()
