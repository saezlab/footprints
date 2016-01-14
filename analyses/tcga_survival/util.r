library(dplyr)
b = import('base')
ar = import('array')
tcga = import('data/tcga')

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

#' A data.frame containing clinical data for all cancer cohorts
clinical = tcga$clinical() %>%
    filter(patient.history_of_neoadjuvant_treatment == "no" &
           is.na(patient.radiations) &
           is.na(patient.follow_ups) &
           is.na(patient.drugs)) %>%
    transmute(study = toupper(admin.disease_code),
              age_days = - as.integer(patient.days_to_birth),
              alive = 1 - as.integer(is.na(patient.days_to_death)),
              surv_days = as.integer(patient.days_to_death %or%
                                     patient.days_to_last_followup),
              barcode = toupper(patient.bcr_patient_barcode),
              sex = as.factor(patient.gender)) %>%
    filter(surv_days > 0) %>% # what is this?
    distinct() %>%
    group_by(barcode) %>%
    filter(age_days == max(age_days))

#' Do pan-cancer survival association using pathway scores and clinical data
#'
#' @param scores  A matrix with samples x pathways
#' @param meta    A data.frame with clinical information
#'                Must have the following fields: surv_days, alive, study, age_days
#' @return        A data.frame with the associations
pancan = function(scores, meta=clinical) {
    ar$intersect(scores, meta$barcode, along=1)
    meta = as.list(meta)
    meta$scores = scores

    #TODO: add sex as covar; but: util tries to subset it, shouldn't
    pancan = st$coxph(surv_days + alive ~ study + age_days + scores,
                      data=meta, min_pts=100) %>%
        filter(term == "scores") %>%
        select(scores, estimate, p.value, size) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))
}

#' Do tissue-specific survival association using pathway scores and clinical data
#'
#' @param scores  A matrix with samples x pathways
#' @param meta    A data.frame with clinical information
#'                Must have the following fields: surv_days, alive, study, age_days
#' @return        A data.frame with the associations
tissue = function(scores, meta=clinical) {
    ar$intersect(scores, meta$barcode, along=1)
    meta = as.list(meta)
    meta$scores = scores

    #TODO: add sex + make it work w/ only one
    tissue = st$coxph(surv_days + alive ~ age_days + scores,
                      subsets=meta$study, data=meta, min_pts=20) %>%
        filter(term == "scores") %>%
        select(scores, subset, estimate, p.value, size) %>%
        group_by(subset) %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        ungroup()
}
