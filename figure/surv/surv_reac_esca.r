library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/reactome.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "surv_reac_esca.pdf"

clinical = tcga$clinical() %>%
    filter(study == "ESCA") %>%
    transmute(age_days = - as.integer(patient.days_to_birth),
              alive = is.na(patient.days_to_death),
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
scores = scores[!nas,grepl("hypoxia", colnames(scores))]
clinical$Hypoxia = "normal"
clinical$Hypoxia[scores > quantile(scores)[4]] = "active"
clinical$Hypoxia[scores < quantile(scores)[1]] = "inactive"
clinical$Hypoxia = factor(clinical$Hypoxia,
    levels = c("normal", "active", "inactive"))
clinical$surv_days = clinical$surv_days / 30.4
clinical = dplyr::filter(clinical, surv_days >= 0)

pdf(OUTFILE, width=10, height=10)
on.exit(dev.off)

library(survival)
library(GGally)
print(survfit(Surv(surv_days, alive) ~ Hypoxia, data=clinical) %>%
    ggsurv() + theme_bw())# + xlim(0, 2000))
