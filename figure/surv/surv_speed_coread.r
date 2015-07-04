library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "surv_speed_coread.pdf"

clinical = tcga$clinical() %>%
    filter(study == "COADREAD") %>%
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
scores = scores[!nas,"H2O2"]
clinical$H2O2 = "normal"
clinical$H2O2[scores > 0.8] = "active"
clinical$H2O2[scores < -0.8] = "inactive"
clinical$H2O2 = factor(clinical$H2O2,
    levels = c("normal", "active", "inactive"))
clinical = dplyr::filter(clinical, surv_days > 1)

pdf(OUTFILE, width=10, height=10)
on.exit(dev.off)

library(survival)
print(survfit(Surv(surv_days, alive) ~ H2O2, data=clinical) %>%
    plt$ggsurv() + xlim(0,1500) + theme_bw())
