library(survival)
library(GGally)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

#' A data.frame containing clinical data for all cancer cohorts
#'
#' Note that we could be filtering by no untreated patients (and tried this),
#' but we think that a treatment that activates or inactivates pathways does
#' impact both pathways and survival and should thus be considered.
clinical = tcga$clinical() %>%
#    filter(patient.history_of_neoadjuvant_treatment == "no" &
#           is.na(patient.radiations) &
#           is.na(patient.follow_ups) &
#           is.na(patient.drugs)) %>%
    transmute(study = toupper(admin.disease_code),
              age_days = - as.integer(patient.days_to_birth),
              alive = 1 - as.integer(is.na(patient.days_to_death)),
              surv_days = as.integer(patient.days_to_death %or%
                                     patient.days_to_last_followup),
              barcode = toupper(patient.bcr_patient_barcode),
              sex = as.factor(patient.gender)) %>%
    mutate(surv_months = surv_days/30.4) %>%
    filter(surv_days > 0) %>% # what is this?
    distinct() %>%
    group_by(barcode) %>%
    filter(age_days == max(age_days)) %>%
    ungroup()

#' Do pan-cancer survival association using pathway scores and clinical data
#'
#' @param scores  A matrix with samples x pathways
#' @param meta    A data.frame with clinical information
#'                Must have the following fields: surv_days, alive, study, age_days
#' @return        A data.frame with the associations
pancan = function(scores, meta=clinical) {
    scores = tcga$map_id(scores, along=1, id_type="patient", subset="primary")
    ar$intersect(scores, meta$barcode, along=1)
    meta = as.list(meta)
    meta$scores = scores

    #TODO: add sex as covar; but: util tries to subset it, shouldn't
    pancan = st$coxph(surv_days + alive ~ study + age_days + scores,
                      data=meta, min_pts=100)

    if (is.data.frame(scores) && all(sapply(scores, is.factor)))
        pancan = pancan %>%
            filter(grepl("^scores", term)) %>%
            mutate(term = sub("scores", "", term)) %>%
            mutate(scores = paste(scores, term, sep="_"),
                   term = "scores")

    pancan %>%
        filter(term == "scores") %>%
        select(scores, estimate, statistic, p.value, size) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))
}

#' Do tissue-specific survival association using pathway scores and clinical data
#'
#' @param scores  A matrix with samples x pathways
#' @param meta    A data.frame with clinical information
#'                Must have the following fields: surv_days, alive, study, age_days
#' @return        A data.frame with the associations
tissue = function(scores, meta=clinical) {
    scores = tcga$map_id(scores, along=1, id_type="patient", subset="primary")
    ar$intersect(scores, meta$barcode, along=1)
    meta = as.list(meta)
    meta$scores = scores

    #TODO: add sex + make it work w/ only one
    tissue = st$coxph(surv_days + alive ~ age_days + scores,
                      subsets=meta$study, data=meta, min_pts=20)

    if (is.data.frame(scores) && all(sapply(scores, is.factor)))
        tissue = tissue %>%
            filter(grepl("^scores", term)) %>%
            mutate(term = sub("scores", "", term)) %>%
            mutate(scores = paste(scores, term, sep="_"),
                   term = "scores")

    tissue %>%
        filter(term == "scores") %>%
        select(scores, subset, estimate, statistic, p.value, size) %>%
        group_by(subset) %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        ungroup()
}

#' Loads a specific TCGA scores file
#'
#' supplies path and extension, and cuts rownames at 16 chars
#'
#' @param id  An identifier, like 'speed_matrix'
#' @return    A matrix with primary tumor (01A) samples as rows and genes and columns
#'            This is mapped to TCGA patient IDs
load_scores = function(id, file=NULL) {
    if (is.null(file))
        file = module_file(io$file_path("../../scores/tcga/pathways_mapped", id, ext=".RData"))
    re =  io$load(file)
    re = re[substr(rownames(re), 14, 16) == "01A",]
    rownames(re) = substr(rownames(re), 1, 12)
    re
}

#' Loads a specific survival association file
#'
#' supplies path and extension, and cuts rownames at 16 chars
#'
#' @param id  An identifier, like 'speed_matrix'
#' @return    A matrix with primary tumor (01A) samples as rows and genes and columns
#'            This is mapped to TCGA patient IDs
load_assocs = function(id, file=NULL) {
    if (is.null(file))
        file = module_file(io$file_path("assocs_cont_mapped", id, ext=".RData"))
    io$load(file)
}

#' Discretizes a numeric vector in "down", "unknown"[x2], "up" characters
#'
#' @param x  A numeric vector
#' @return   A character vector
discretize_quartiles = function(x, numeric=TRUE, na.rm=TRUE) {
    qq = quantile(x, na.rm=na.rm)
    re = rep(0, length(x))
    re[x > unname(qq[4])] = 1
    re[x < unname(qq[2])] = -1

    if (!numeric) {
        re = as.character(re)
        re[re == "-1"] = "down"
        re[re == "0"] = "unknown"
        re[re == "1"] = "up"
    }
    re
}

#' Converts a row of the associations data.frame to a survival fit
#'
#' @param row     A data.frame row with the fields 'subset' and 'adj.p'
#' @param rename  A named vector which levels should have which names
#'                The names must be (-1,0,1) and the values the names
row2survFit = function(row, rename=c("-1"="down", "0"="unknown", "1"="up")) {
    score = scores
    clin = clinical
    if ("subset" %in% names(row))
        clin = filter(clinical, study == row['subset'])

    ar$intersect(score, clin$barcode, along=1)
    clin$pathway = score[,sub("_.*$", "", row['scores'])]

    if (!is.null(rename))
        clin$pathway = rename[as.character(clin$pathway)]

    survfit(Surv(surv_months, alive) ~ pathway, data=clin) %>%
        ggsurv() +
            xlim(0, 52) +
            theme_bw() +
            xlab("Survival (weeks)") +
            ggtitle(paste(row['subset'], row['scores'], row['adj.p']))
}
