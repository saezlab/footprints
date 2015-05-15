library(dplyr)
b = import('base')
io = import('io')
st = import('stats')
plt = import('plot')
icgc = import('icgc')

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

scores = io$load(OUTFILE)
survival = clinical[c('known_survival_time','donor_vital_status')]
rownames(survival) = clinical$icgc_sample_id
#TODO: make sure they are ordered the same

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

# run cox regression for associations (tissue as covariate)
st$cox(survival ~ scores, subsets=clinical$tissue) %>%
    mutate(p.adj = p.adjust(p.value, method="fdr") %>%
    plt$volcano() %>%
    print()
