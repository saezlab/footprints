io = import('io')
roc = import('./roc')

sdata = io$load('sigs_scores.RData')

debug(roc$scores2roc)
lookup = setNames(sub("\\..*$", "", colnames(sdata$scores)),
                  colnames(sdata$scores))
sroc = roc$scores2roc(sdata, lookup)
