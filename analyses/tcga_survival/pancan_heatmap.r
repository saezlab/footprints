library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')

OUTFILE = commandArgs(TRUE)[2] %or% "pancan_heatmap.pdf"

clinical = tcga$clinical() %>%
    transmute(study = study,
              age_days = - as.integer(patient.days_to_birth),
              alive = 1 - as.integer(is.na(patient.days_to_death)),
              surv_days = as.integer(patient.days_to_death %or%
                                     patient.days_to_last_followup),
              barcode = toupper(patient.bcr_patient_barcode),
              gender = as.factor(patient.gender)) %>%
    filter(surv_days > 0) #FIXME:

file_to_assocs = function(fname, clinical) {
    scores = io$load(fname)
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
    clinical = as.list(clinical)
    clinical$scores = scores

    # tissue covariate
    st$coxph(surv_days + alive ~ study + age_days + scores, 
             data=clinical, min_pts=100) %>%
        filter(term == "scores") %>%
        select(scores, estimate, p.value, size) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))
}

fnames = c("speed_matrix", "gsea_reactome", "gsea_go", "spia", "pathifier")
assocs = lapply(fnames, function(fn) {
    fname = paste0("../../scores/tcga/", fn, ".RData")
    file_to_assocs(fname, clinical) %>%
        mutate(method = fn)
}) %>%
    bind_rows() %>%
    mutate(estimate = ifelse(adj.p < 0.1, estimate, NA)) %>%
    mutate(label = ifelse(adj.p < 0.05, "*", "")) %>%
    mutate(label = ifelse(adj.p < 0.01, "***", label))

pdf(OUTFILE, paper="a4r", width=10, height=5)
plt$matrix(assocs, estimate ~ method + scores) +
    xlab("Pathway") +
    ylab("Method") +
    ggtitle("Pan-cancer survival associations")
dev.off()
