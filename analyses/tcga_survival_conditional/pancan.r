library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_norm.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "cont_speed_norm.pdf"

clinical = tcga$clinical() %>%
    transmute(study = study,
              age_days = - as.integer(patient.days_to_birth),
              alive = 1 - as.integer(is.na(patient.days_to_death)),
              surv_days = as.integer(patient.days_to_death %or%
                                     patient.days_to_last_followup),
              barcode = toupper(patient.bcr_patient_barcode),
              gender = as.factor(patient.gender)) %>%
    filter(surv_days > 0) #FIXME:

scores = io$load(INFILE)
# subset primary tumors only
#TODO: what if i do average over all tumor samples for each patient
#  (but strong signals should be in there either way)
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

if (nrow(clinical) < 10)
    stop("survival+expression scores < 10 observations, stopping")

clinical = as.list(clinical)
clinical$scores = scores

pdf(OUTFILE, paper="a4r", width=26, height=20)

# tissue covariate
#TODO: add gender as covar; but: util tries to subset it, shouldn't
st$coxph(surv_days + alive ~ study + age_days + scores, data=clinical, min_pts=100) %>%
    filter(term == "scores") %>%
    select(scores, estimate, p.value, size) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1) %>%
    print()

# separate regressions for each tissue
#TODO: add gender + make it work w/ only one
st$coxph(surv_days + alive ~ age_days + scores, subsets=clinical$study, data=clinical, min_pts=20) %>%
    filter(term == "scores") %>%
    select(scores, subset, estimate, p.value, size) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup() %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) + ggtitle(sum(clinical$adj.p < 0.1)) %>%
    print()

dev.off()
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
st = import('stats')
gdsc = import('data/gdsc')
plt = import('plot')

scores2assocs = function(scores, fdr=0.1) {
    assocs = st$lm(Ys ~ scores, subsets=tissues) %>%
        filter(term == "scores") %>%
        group_by(subset) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        ungroup() %>%
        mutate(label = paste(Ys, scores, subset, sep=":"))
    message(sum(assocs$adj.p < fdr), " significant associations, FDR: ", fdr)
    assocs
}

conditional_assocs = function(scores, fdr=0.1) {
    best_speed_hit = a_speed %>%
        group_by(Ys, subset) %>%
        filter(adj.p == min(adj.p)) %>%
        select(Ys, scores, subset, adj.p) %>%
        distinct() # duplicate drug names in Ys

    condition = ar$like(NA, Ys)
    for (i in 1:nrow(best_speed_hit)) {
        cur = best_speed_hit[i,]
        condition[tissues == cur$subset, cur$Ys] =
            s_speed[tissues == cur$subset, cur$scores]
    }

    # return df with associations still significant after conditioning
    st$lm(Ys ~ condition + scores, group=c("Ys", "condition"), subsets=tissues) %>%
        filter(term == "scores") %>%
        group_by(subset) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        ungroup() %>%
        filter(adj.p < fdr)
}

assocs2plot = function(assocs, ylim, fdr=0.1, condition=NULL) {
    assocs = assocs %>%
         plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1, thresh=fdr)

    if (!is.null(condition)) {
        subset_valid = condition %>%
            select(Ys, scores) %>%
            df$contains(assocs, .) %catch% FALSE

        assocs[!subset_valid & assocs$adj.p < fdr,'color'] = "#00000000"
        assocs[!subset_valid & assocs$adj.p < fdr,'circle'] = "#767676ff"
    }

    assocs %>%
        plt$volcano(base.size=2, p=fdr, xlim=c(-5,5), ylim=ylim)
}

# load scores; s_: prefix for scores
s_speed = io$load('../../scores/gdsc/speed_matrix.RData')
s_reac = io$load('../../scores/gdsc/reactome.RData')
s_go = io$load('../../scores/gdsc/go.RData')

# load sanger data
Ys = gdsc$drug_response('IC50s', min_tissue_measured=2) # or AUC
tissues = gdsc$tissues(minN=5)
ar$intersect(s_speed, s_reac, s_go, tissues, Ys, along=1)

# a_: prefix for associations
a_speed = scores2assocs(s_speed, 0.1) # 0.2: 190, 0.1: 74
a_reac = scores2assocs(s_reac, 0.1) # 0.2: 130, 0.1: 36
a_go = scores2assocs(s_go, 0.1) # 0.2: 66, 0.1: 12

# c_: prefix conditional associations
c_reac = conditional_assocs(s_reac)
c_go = conditional_assocs(s_go)

# v_: prefix for volcano plots
v_speed = assocs2plot(a_speed, c(1,1e-6)) # limit are there so axes are equal
v_reac = assocs2plot(a_reac, c(1,1e-6), condition=c_reac)
v_go = assocs2plot(a_go, c(1,1e-6), condition=c_go)

# arrange plots
pcol = arrangeGrob(v_speed, v_reac, v_go, ncol=1, nrow=3)

# save to pdf
pdf("tissue.pdf", width=6, height=15)
plot(pcol)
dev.off()
