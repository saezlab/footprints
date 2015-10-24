library(ggplot2)
library(dplyr)
b = import('base')
io = import('io')

null = io$load("model_vs_mut_NULL_100rep_alldrugs.RData") %>%
    transmute(drugs = drugs,
              dset = dset,
              subset = subset,
              rep = rep,
              mse.test.null = mse.test.mean) %>%
    mutate(mse.test.null[is.na(mse.test.null)] = Inf)

model = io$load("model_vs_mut_alldrugs.RData") %>%
    select(-rmse.test.sqrt.of.mean, -mae.test.mean) %>%
    filter(drugs %in% null$drugs) %>%
    left_join(null) %>%
    group_by(drugs, dset, subset) %>%
    summarize(emp.p = (1 + sum(mse.test.mean > mse.test.null)) / n(),
              mse.test.mean = mse.test.mean[1]) %>%
    mutate(adj.p = p.adjust(emp.p, method="fdr"))

model %>%
    filter(adj.p < 0.1) %>%
    ggplot(aes(dset)) + geom_bar()
