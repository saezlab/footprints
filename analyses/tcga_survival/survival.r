library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
icgc = import('data/icgc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

clinical = icgc$clinical() %>% #TODO: select tumor only
#    filter(tumour_confirmed == "yes") %>%
#    filter(specimen_donor_treatment_type == "no treatment") %>%
    mutate(known_survival_time = donor_survival_time %or% donor_interval_of_last_followup) %>%
    select(icgc_specimen_id, known_survival_time, donor_vital_status,
           donor_sex, donor_age_at_diagnosis, tissue) %>%
    filter(!is.na(known_survival_time),
           donor_vital_status %in% c('alive','deceased'))

# possible questions here:
#  using all tumor data, is pathway activity associated with survival outcome?
#    - subset treatment naive?
#    - can it predict relapse?
#    - does a treatment activate pathways?
# -- all in covariate and subset tissue data

scores = io$load(INFILE)
rownames(clinical) = clinical$icgc_specimen_id
ar$intersect(scores, clinical, along=1) #TODO: check if works with data.frame; if not, make it work

clinical = as.list(clinical)
clinical$donor_vital_status = as.integer(clinical$donor_vital_status == "alive")
clinical$scores = scores

if (nrow(clinical) < 10)
    stop("survival+expression scores < 10 observations, stopping")

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off())

# tissue covariate
st$coxph(time + status ~ donor_sex + tissue + scores) %>%
    filter(term == "scores") %>%
    select(-time, -status, -tissue, -term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    plt$color$p_effect() %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1) %>%
    print()

# separate regressions for each tissue
st$coxph(time + status ~ donor_sex + scores, subsets=tissue) %>%
    select(-time, -status, -term) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup() %>%
    plt$color$p_effect() %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) %>%
    print()
