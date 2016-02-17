# point of this file:
# - use the zscores to create a linear model
library(dplyr)
b = import('base', attach_operators=FALSE)
import('base/operators')
io = import('io')
ar = import('array')
st = import('stats')

#' Fits a linear model on Z-scores
#'
#' @param zdata  A list with the zscore matrix and index object
#' @return       The coefficients matrix [gene x pathway]
zscore2model = function(zdata, hpc_args=NULL) {
	index = zdata$index
	zscores = t(zdata$zscores) * index$sign

	# fit model to pathway perturbations
	pathway = t(ar$mask(index$pathway)) + 0
	pathway["EGFR",] = pathway["EGFR",] + pathway["MAPK",] + pathway["PI3K",]
	pathway["TNFa",] = pathway["TNFa",] + pathway["NFkB",]

	mod = st$lm(zscores ~ 0 + pathway, data=index, min_pts=30, atomic="pathway",
				hpc_args=hpc_args) %>%
		transmute(gene = zscores,
				  pathway = sub("^pathway", "", term),
				  zscore = estimate,
				  p.value = p.value) %>%
		group_by(gene) %>%
		mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
		ungroup()

	zfit = ar$construct(zscore ~ gene + pathway, data=mod)
	pval = ar$construct(p.adj ~ gene + pathway, data=mod)

    zfit
}

## best x-val: 0.15, 1e-4 -> use which?
#zfit[abs(zfit) < 1 | pval > 1e-3 | is.na(zfit)] = 0
#print(colSums(zfit != 0))

if (is.null(module_name())) {
    ZDATA = commandArgs(TRUE)[1] %or% "../data/zscores.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "model_linear.RData"

    # load speed data, index; filter for train set only
    zdata = io$load(ZDATA)
    zfit = zscore2model(zdata, hpc_args=list(n_jobs=10, memory=2048))

    # save resulting object
    save(zfit, file=OUTFILE)
}
