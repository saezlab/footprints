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
	mod = st$lm(zscores ~ 0 + pathway, data=index, min_pts=100,
				hpc_args=hpc_args) %>%
		transmute(gene = zscores,
				  pathway = sub("^pathway", "", term),
				  zscore = estimate,
				  p.value = p.value) %>%
		mutate(adj.p = p.adjust(p.value, method="fdr"))

	zfit = ar$construct(zscore ~ gene + pathway, data=mod)
	pval = ar$construct(p.value ~ gene + pathway, data=mod)

	# filter zfit to only include top 100 genes per pathway
    model = zfit
	model[apply(pval, 2, function(p) !b$min_mask(p, 100))] = 0

    list(assocs=mod, model=model)
}

if (is.null(module_name())) {
	ZDATA = commandArgs(TRUE)[1] %or% "../data/zscores.RData"
	OUTFILE = commandArgs(TRUE)[2] %or% "model_linear.RData"

	# load speed data, index; filter for train set only
	zdata = io$load(ZDATA)
	result = zscore2model(zdata, hpc_args=list(n_jobs=10, memory=2048))

	# save resulting object
	save(result, file=OUTFILE)
}
