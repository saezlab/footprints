# point of this file:
# - use the dscores to create a linear model
library(dplyr)
b = import('base', attach_operators=FALSE)
import('base/operators')
io = import('io')
ar = import('array')
st = import('stats')

INFILE = commandArgs(TRUE)[1] %or% "../data/scores.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "model_ridge.RData"

# load speed data, index; filter for train set only
zobj = io$load('../data/scores.RData')
index = zobj$index #%>% filter(is.na(exclusion))
dscores = zobj$dscores[,index$id]

# adjust object for linear modelling
inh = index$effect=="inhibiting"
dscores[,inh] = -dscores[,inh]
dscores = t(dscores)

do_path = function(path) {
    path_sign = ifelse(index$effect == "activating", 1, -1)
    path_mask = index$pathway == path
    if (path == "NFkB")
        path_mask = index$pathway %in% c("NFkB", "TNFa")
    if (path == "MAPK")
        path_mask = index$pathway %in% c("MAPK", "EGFR")
    if (path == "PI3K")
        path_mask = index$pathway %in% c("PI3K", "EGFR")
    path_sign[!path_mask] = 0

#    re = st$ml(path_sign ~ dscores, train_args=list("regr.glmnet", dfmax=100, alpha=0.5),
#               xval=5, models=TRUE, result_only=TRUE)[[1]]

    library(glmnet)
    re = cv.glmnet(x=dscores, y=path_sign, nfolds=5, alpha=0.05, dfmax=100)
    coeff = coef(re, s="lambda.1se")
    coeff = setNames(coeff[,1], rownames(coeff))[-1]
}

# fit model to pathway perturbations
zfit = sapply(unique(index$pathway), do_path, USE.NAMES=TRUE, simplify=FALSE) %>%
    ar$stack(along=2)
zfit = zfit[,colSums(zfit!=0) != 0]

# save resulting object
save(zfit, file=OUTFILE)
