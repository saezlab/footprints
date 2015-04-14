library(dplyr)
b = import('base')
icgc = import('icgc')

# load clinical data, extract survival time+status
clinicalsample = icgc$getClinicalSample() %>%
    select(icgc_specimen_id, icgc_sample_id)

clinical = icgc$getClinical() %>% #TODO: select tumor only
#    filter(tumour_confirmed == "yes") %>%
#    filter(specimen_donor_treatment_type == "no treatment") %>%
    mutate(known_survival_time = donor_survival_time %or% donor_interval_of_last_followup) %>%
    select(icgc_specimen_id, known_survival_time, donor_vital_status,
           donor_sex, donor_age_at_diagnosis, tissue) %>%
    filter(!is.na(known_survival_time),
           donor_vital_status %in% c('alive','deceased')) %>%
    left_join(clinicalsample, by="icgc_specimen_id") %>%
    filter(icgc_sample_id %in% icgc$namesRNASeq()[[1]])

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - filter for last known alive? [where is this field?]
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

save(clinical, file="survival.RData")
