library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
hpc = import('hpc')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/pathways_mapped/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_matrix.RData"

#' Load required data
data = list(
    scores = io$load(INFILE),
    response = gdsc$drug_response('IC50s'),
    clinical = gdsc$drug_response('IC50s', min_tissue_measured=0, stage=2),
    noexp = gdsc$drug_response('IC50s', min_tissue_measured=0, median_top=10, stage=1),
    sensi = gdsc$drug_response('IC50s', min_tissue_measured=5, median_top=10),
    clin_sens = gdsc$drug_response('IC50s', min_tissue_measured=0, stage=2, median_top=10),
    tissues = gdsc$tissues(),
    MSI = (gdsc$MASTER_LIST['MMR'] == "MSI-H")[,1]
) %>% ar$intersect_list(along=1)

#' Tissues as covariate
#'
#' @param data      A list with: tissue, MSI, scores
#' @return          A data.frame with the association results
assocs.pan = st$lm(response ~ tissues + MSI + scores, data=data, min_pts=50,
                   hpc_args = list(job_size = 5000)) %>%
    filter(term == "scores") %>%
    mutate(subset = tissue) %>%
    select(-term, subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

#' Tissues as subsets
#'
#' @param resp_sub  A character string of which `response` to subset
#' @param data      A list with: tissue, MSI, scores
#' @return          A data.frame with the association results
tissue_assocs = function(resp_sub, data) {
    st = import('stats')
    data$response = data[[resp_sub]]
    re = st$lm(response ~ MSI + scores, subsets=data$tissues, data=data, min_pts=10,
               hpc_args = list(job_size = 5000)) %>%
        mutate(tissue = subset, subset=resp_sub)
}

assocs.tissue = sapply(c('clinical', 'noexp', 'sensi'), tissue_assocs, data=data, simplify=FALSE) %>%
    bind_rows() %>%
    filter(term == "scores") %>%
    select(-term) %>%
    group_by(tissue, subset) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

save(assocs.pan, assocs.tissue, file=OUTFILE)
