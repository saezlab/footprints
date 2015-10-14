library(dplyr)
library(reshape2)
library(ggplot2)
io = import('io')
gdsc = import('data/gdsc')
st = import('stats')

tissues = gdsc$tissues(minN=10)

muts = gdsc$mutated_genes()[,c("HRAS", "KRAS", "MAP2K1", "BRAF")] %>%
    melt() %>%
    filter(value == TRUE) %>%
    select(-value) %>%
    transmute(cosmic = as.character(Var1), gene = Var2)

speed = io$load('../../scores/gdsc/speed_matrix.RData')[,"MAPK"] %>%
    data.frame(cosmic=names(.),
               mapk = .,
               up = . > quantile(.)[4],
               down = . < quantile(.)[2]) %>%
    select(-mapk) %>%
    melt(id="cosmic") %>%
    filter(value == TRUE) %>%
    transmute(cosmic=cosmic, mapk=variable)

MEKis = c("AZ628", "Dabrafenib", "(5Z)-7-Oxozeaenol",
          "PD-0325901", "RDEA119", "Trametinib", "CI-1040", "VX-11e")
resp = gdsc$drug_response()[,MEKis] %>%
    melt() %>%
    transmute(cosmic = as.character(Var1), drug = Var2, IC50 = value) %>%
    left_join(muts) %>%
    left_join(speed) %>%
    mutate(no_mut = is.na(gene))
resp$no_mapk_up = resp$no_mut | resp$mapk!="up" | is.na(resp$mapk)
resp$no_both = !resp$no_mut | !resp$no_mapk_up

pdf("pancan.pdf", width=10, height=8)

ggplot(resp, aes(IC50, fill=no_mut)) +
    stat_density(aes(y = ..count..), position="stack", color="black") +
    theme_bw()

ggplot(resp, aes(IC50, fill=no_mapk_up)) +
    stat_density(aes(y = ..count..), position="stack", color="black") +
    theme_bw()

ggplot(resp, aes(IC50, fill=no_both)) +
    stat_density(aes(y = ..count..), position="identity", color="black") +
    theme_bw()

dev.off()
