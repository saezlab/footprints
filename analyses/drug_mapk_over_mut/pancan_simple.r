library(dplyr)
library(reshape2)
library(ggplot2)
io = import('io')
gdsc = import('data/gdsc')
st = import('stats')

tissues = gdsc$tissues(minN=10)

drug2plot = function(drug, pathway="MAPK", genes = c("HRAS", "KRAS", "BRAF")) {
    muts = gdsc$mutated_genes()[,genes,drop=FALSE] %>%
        apply(1, any)
    muts = data_frame(
        cosmic = as.character(names(muts)),
        mapk = ifelse(muts, "mut", "wt")
    )

    if (pathway == "p53")
        muts$mapk = ifelse(muts$mapk=="mut","wt","mut")

    speed = io$load('../../scores/gdsc/speed_matrix.RData')[,pathway] %>%
        data_frame(cosmic=names(.),
                   mapk = .,
                   up = . > quantile(.)[4],
                   down = . < quantile(.)[2]) %>%
#                   average = !(up | down)) %>% # alternatively: 4,2 <-> 3,3
        select(-mapk) %>%
        melt(id="cosmic") %>%
        filter(value == TRUE) %>%
        transmute(cosmic=cosmic, mapk=as.character(variable))

    both = muts %>%
        mutate(gene = mapk) %>%
        select(-mapk) %>%
        inner_join(speed) %>%
        transmute(cosmic=cosmic, mapk=ifelse(gene=="mut"&mapk=="up", "both", "other"))

    combined = bind_rows(muts, speed, both)

    resp = gdsc$drug_response()[,drug] %>%
        as.data.frame() %>%
        add_rownames() %>%
        transmute_(cosmic = as.character("rowname"), IC50=".") %>%
        right_join(combined)

    resp$mapk = factor(resp$mapk, levels=c("wt", "mut", "down", "average", "up", "both"))
    resp = na.omit(resp)

    ggplot(resp, aes(x=mapk, y=IC50)) +
        geom_boxplot(outlier.shape = NA) +
        theme_bw() +
        ggtitle(paste0(drug, " [", paste(genes, collapse=", "), "]"))
}

MEKis = c("AZ628", "Dabrafenib", "(5Z)-7-Oxozeaenol",
      "PD-0325901", "RDEA119", "Trametinib", "CI-1040", "VX-11e")

pdf("pancan.pdf", width=10, height=8)

dev.off()
