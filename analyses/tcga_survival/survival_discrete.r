library(survival)
library(GGally)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

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

scores = io$load(INFILE)
# subset primary tumors only
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

scores = ar$map(scores, along=1, subsets=clinical$study, function(x) {
    qq = quantile(x)
    re = rep(0, length(x))
    re[x > unname(qq[4])] = 1
    re[x < unname(qq[2])] = -1
    if (sum(abs(re)) < 10)
        rep(NA, length(x))
    else
        re
})

if (nrow(clinical) < 10)
    stop("survival+expression scores < 10 observations, stopping")

# tissue covariate
#assocs.pan = st$coxph(surv_days + alive ~ study + gender + age_days + scores,
assocs.pan = st$coxph(surv_days + alive ~ study + scores,
        data=clinical, min_pts=100) %>%
    filter(term == "scores") %>%
    select(scores, estimate, p.value, size) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

# separate regressions for each tissue
#assocs.tissue = st$coxph(surv_days + alive ~ gender + age_days + scores,
assocs.tissue = st$coxph(surv_days + alive ~ scores,
        subsets=clinical$study, data=clinical, min_pts=20) %>%
    filter(term == "scores") %>%
    select(scores, subset, estimate, p.value, size) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

row2survFit = function(row, include_normal=FALSE) {
    if ("subset" %in% names(row))
        study_filter = clinical$study == row['subset']
    else
        study_filter = rep(TRUE, nrow(clinical))

    clin = clinical[study_filter,]
    pathway = as.factor(scores[study_filter, row['scores']])
    stopifnot(levels(pathway) == c("-1", "0", "1"))
    levels(pathway) = c("inactive", "normal", "active")

    if (!include_normal)
        pathway[pathway == "normal"] = NA

    clin$surv_months = clin$surv_days / 30.4
    clin$pathway = pathway

    survfit(Surv(surv_months, alive) ~ pathway, data=clin) %>%
        ggsurv() +
            xlim(0, 52) +
            theme_bw() +
            xlab("Survival (weeks)") +
            ggtitle(paste(row['subset'], row['scores'], row['adj.p'])) %>%
        print()
}

pdf(OUTFILE, paper="a4r", width=26, height=20)

assocs.pan %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1) %>%
    print()

assocs.pan %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(5) %>%
    apply(1, row2survFit)

assocs.tissue %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) %>%
    print()

assocs.tissue %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(15) %>%
    apply(1, row2survFit)

dev.off()
