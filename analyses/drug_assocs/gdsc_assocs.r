library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/pathways_mapped/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.RData"

# load required data
scores = io$load(INFILE)
Ys = gdsc$drug_response('IC50s') # or AUC
Ys_clinical = gdsc$drug_response('IC50s', min_tissue_measured=0, stage=2)
Ys_noexp = gdsc$drug_response('IC50s', min_tissue_measured=0, median_top=10, stage=1)
Ys_sensi = gdsc$drug_response('IC50s', min_tissue_measured=5, median_top=10)
Ys_clin_sens = gdsc$drug_response('IC50s', min_tissue_measured=0, stage=2, median_top=10)
tissues = gdsc$tissues(minN=15)
ar$intersect(scores, tissues, Ys, Ys_clinical, Ys_noexp, Ys_sensi, Ys_clin_sens, along=1)

# tissues as covariate
assocs.pan = st$lm(Ys ~ tissues + scores, hpc_args=list(job_size=2000)) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

# tissues as subsets
assocs.tissue = function(Yf) st$lm(Yf ~ scores, subsets=tissues, hpc_args=list(job_size=10000)) %>%
    filter(term == "scores") %>%
    select(-term) %>%
    group_by(subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

assocs.tissue_clinical = assocs.tissue(Ys_clinical)
assocs.tissue_noexp = assocs.tissue(Ys_noexp)
assocs.tissue_sensi = assocs.tissue(Ys_sensi)
assocs.tissue_clin_sens = assocs.tissue(Ys_clin_sens)

save(assocs.pan, assocs.tissue_clinical, assocs.tissue_noexp, assocs.tissue_sensi, assocs.tissue_clin_sens, file=OUTFILE)
