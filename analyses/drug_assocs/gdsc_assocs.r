library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
hpc = import('hpc')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/pathways_mapped/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.RData"

# load required data
data = list(
    scores = io$load(INFILE),
    Ys = gdsc$drug_response('IC50s'),
    clinical = gdsc$drug_response('IC50s', min_tissue_measured=0, stage=2),
    noexp = gdsc$drug_response('IC50s', min_tissue_measured=0, median_top=10, stage=1),
    sensi = gdsc$drug_response('IC50s', min_tissue_measured=5, median_top=10),
    clin_sens = gdsc$drug_response('IC50s', min_tissue_measured=0, stage=2, median_top=10),
    tissues = gdsc$tissues(minN=15)
) %>% ar$intersect_list(along=1)

scores = data$scores
tissues = data$tissues
Ys = data$Ys
Ysubs = data[c('clinical', 'noexp', 'sensi', 'clin_sens')]

#' Tissues as covariate
#'
#' @param scores    A matrix of [cell lines x set or pathway scores]
#' @param Ys        A matrix drug response [cell lines x drugs]
#' @param scoresub  A character string of which column of `scores` to subset
#' @param tissues   A character vector of tissues, same cell lines as `scores`
#' @return          A data.frame with the association results
pan_assocs = function(scores, Ys, scoresub, tissues) {
    ar = import('array')
    st = import('stats')

    score = scores[, scoresub, drop=FALSE]

    assocs.pan = st$lm(Ys ~ tissues + score, data=data) %>%
        filter(term == "score") %>%
        select(-term)
}

#' Tissues as subsets
#'
#' @param scores    A matrix of [cell lines x set or pathway scores]
#' @param Ysubs     A list of different subsets for drug response
#' @param scoresub  A character string of which column of `scores` to subset
#' @param ysub      A character string of which `Ysubs` to subset
#' @param tissues   A character vector of tissues, same cell lines as `scores`
#' @return          A data.frame with the association results
tissue_assocs = function(scores, Ysubs, scoresub, ysub, tissues) {
    ar = import('array')
    st = import('stats')

    data = list(Yf = Ysubs[[ysub]], score = scores[,scoresub,drop=FALSE]) %>%
        ar$intersect_list(along=1)

    re = st$lm(Yf ~ score, subsets=tissues, data=data) %>%
        filter(term == "score") %>%
        select(-term) %>%
        mutate(Ysub = ysub)
}

assocs.pan = hpc$Q(pan_assocs, scoresub=colnames(scores),
                   const = list(scores=scores, Ys=Ys, tissues=tissues),
                   job_size = 100) %>%
    bind_rows() %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

assocs.tissue = hpc$Q(tissue_assocs, scoresub=colnames(scores), ysub=names(Ysubs),
                      const = list(scores=scores, Ysubs=Ysubs, tissues=tissues),
                      expand_grid=TRUE, job_size=100) %>%
    bind_rows() %>%
    group_by(subset, Ysub) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

save(assocs.pan, assocs.tissue, file=OUTFILE)
