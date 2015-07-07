library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/reactome.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "surv_reac_luad.pdf"

clinical = tcga$clinical() %>%
    filter(study == "LUAD") %>%
    transmute(age_days = - as.integer(patient.days_to_birth),
              alive = 1 - as.integer(is.na(patient.days_to_death)),
              surv_days = as.integer(patient.days_to_death %or%
                                     patient.days_to_last_followup),
              barcode = toupper(patient.bcr_patient_barcode),
              gender = as.factor(patient.gender))

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

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
nas = is.na(cc$alive)
clinical = cc[!nas,]
scores = scores[!nas,grepl("Wnt", colnames(scores))]
clinical$Wnt = "normal"
clinical$Wnt[scores > quantile(scores)[4]] = "active"
clinical$Wnt[scores < quantile(scores)[1]] = "inactive"
clinical$Wnt = factor(clinical$Wnt,
    levels = c("normal", "active", "inactive"))
clinical$surv_days = clinical$surv_days / 30.4
clinical = dplyr::filter(clinical, surv_days >= 0)

pdf(OUTFILE, width=8, height=6)
on.exit(dev.off)

library(survival)
print(survfit(Surv(surv_days, alive) ~ Wnt, data=clinical) %>%
    plt$ggsurv() + theme_bw())# + xlim(0,2000))
