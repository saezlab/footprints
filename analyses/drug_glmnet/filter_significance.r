#library(fitdistrplus)
library(dplyr)
b = import('base')
io = import('io')

#null = io$load("model_vs_mut_NULL_1000rep_20drugs.RData")
null = io$load("model_vs_mut_NULL_100rep_alldrugs.RData") %>%
    group_by(drugs, dset, subset) %>%
    do(data.frame(mse.cutoff = b$minN(.$mse.test.mean, 5)))

model = io$load("model_vs_mut_alldrugs.RData") %>%
    filter(drugs %in% null$drugs) %>%
    inner_join(null) %>%
    mutate(valid = mse.test.mean < mse.cutoff) %>%
    group_by(drugs, subset) %>%
    filter(valid) %>%
    mutate(best_model = mse.test.mean == min(mse.test.mean))

# total number of models
# significant number of models
# plot: # of best significant models

#drug = sample(model$drugs, 1)
#method = "speed"
#tissue = "BRCA"
#nulldist = null %>%
#    filter(drugs == drug & dset == method & subset == tissue)
#
# this does not really follow any distribution
# can do it empirically, will need a lot of null models
#dd = nulldist$mse.test.mean
#hist(dd, 100)
#plot(fitdist(dd[!is.na(dd)], "norm"))

# select all models that have emp_p<0.05
# can i just group+merge?
