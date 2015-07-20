library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
tcga = import('data/tcga')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "../../expr_cluster/speed_cluster.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "surv_speed.pdf"

# load scores, from tumors only
clusters = io$load(INFILE)
clusters = clusters[lapply(clusters, length) != 0]
clusters$SKCM = NULL #FIXME: only cell lines in there
clusters$THCA = NULL #FIXME: only cell lines in there
scores = lapply(clusters, function(cl) {
    re = cl %>%
        lapply(function(x) x[grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-0[13]A$", rownames(x)),]) %>%
        ar$stack(along=2)
    rownames(re) = substr(rownames(re), 1, 12)
    re[,colSums(re)>1]
})

# save pdf w/ pan-cancer & tissue specific
pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

#TODO: this is quite similar code to analyses/tcga_survival_surv_assocs
for (tissue in names(scores)) {
    message(tissue)
    score = scores[[tissue]]
    clinical = tcga$clinical(tissue = tissue) %>%
        transmute(age_days = - as.integer(patient.days_to_birth),
                  alive = 1 - as.integer(is.na(patient.days_to_death)),
                  surv_days = as.integer(patient.days_to_death %or%
                                         patient.days_to_last_followup),
                  barcode = toupper(patient.bcr_patient_barcode),
                  gender = as.factor(patient.gender)) %>%
        distinct()

    rownames(clinical) = clinical$barcode
    ar$intersect(score, clinical)

    if (nrow(clinical) == 0)
        next

    #TODO: add gender if both present
    re = st$coxph(surv_days + alive ~ age_days + score, data=clinical, min_pts=20) %>%
        filter(term == "score") %>%
        select(score, estimate, p.value, size) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               label = score) %>%
        plt$color$p_effect("adj.p") %>%
        plt$volcano(p=0.1) + ggtitle(tissue)
    print(re)
}
