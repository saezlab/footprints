library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/pathways_mapped/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "assocs_mapped/speed_matrix.RData"

#' Load required data
data = list(
    scores = io$load(INFILE),
    drug = gdsc$drug_response('IC50s'),
    clinical = gdsc$drug_response('IC50s', min_tissue_measured=0, stage=2),
    noexp = gdsc$drug_response('IC50s', min_tissue_measured=0, median_top=10, stage=1),
    sensi = gdsc$drug_response('IC50s', min_tissue_measured=5, median_top=10),
#    clin_sens = gdsc$drug_response('IC50s', min_tissue_measured=0, stage=2, median_top=10),
    tissues = gdsc$tissues(),
    MSI = (gdsc$MASTER_LIST['MMR'] == "MSI-H")[,1]
) %>% ar$intersect_list(along=1)

#' Tissues as covariate
#'
#' @param data      A list with: tissue, MSI, scores
#' @return          A data.frame with the association results
pan = st$lm(drug ~ tissues + MSI + scores, data=data, min_pts=50,
            hpc_args = list(job_size=1e4, n_jobs=200, memory=20480)) %>%
    tbl_df() %>%
    filter(term == "scores") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

#' Tissues as subsets
#'
#' @param resp_sub  A character string of which `drug` to subset
#' @param data      A list with: tissue, MSI, scores
#' @return          A data.frame with the association results
tissue_assocs = function(resp_sub="all", data=get("data", envir=parent.env(environment()))) {
    st = import('stats')
    data$drug = data[[resp_sub]]
    gc() # something is not properly cleaned up here
    re = st$lm(drug ~ MSI + scores, subsets=data$tissues, data=data, min_pts=10,
               hpc_args = list(job_size=1e4, n_jobs=300, memory=20480)) %>%
        mutate(tissue = subset, subset=resp_sub)
    gc() # something is not properly cleaned up here

    # we are only interested in pathway scores here, not MSI
    re %>%
        filter(term == "scores") %>%
        select(-term) %>%
        group_by(tissue, subset) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        ungroup()
}

#' Tissues as subsets
#'
#' sensi @GO: long vectors not supported w/ send_common_data
#'
#' @return  A list of data frames for subsets clinical, noexp and sensi
tissue = c('all', 'clinical', 'noexp', 'sensi') %>%
    b$lnapply(function(subs) tissue_assocs(subs, data=data))

save(pan, tissue, file=OUTFILE)
