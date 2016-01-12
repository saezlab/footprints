library(survival)
library(GGally)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')
tcga = import('data/tcga')
util = import('./util')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

discretize_quartiles = function(x) {
    qq = quantile(x)
    re = rep(0, length(x))
    re[x > unname(qq[4])] = 1
    re[x < unname(qq[2])] = -1
    if (sum(abs(re)) < 10)
        rep(NA, length(x))
    else
        re
}

row2survFit = function(row, include_normal=FALSE) {
    if ("subset" %in% names(row))
        study_filter = clinical$study == row['subset']
    else
        study_filter = rep(TRUE, nrow(clinical))

    clin = clinical[study_filter,]
    pathway = as.factor(scores[study_filter, row['scores']])
    stopifnot(levels(pathway) == c("-1", "0", "1"))
    levels(pathway) = c("inactive", "normal", "active")

    if (!include_normal)
        pathway[pathway == "normal"] = NA

    clin$surv_months = clin$surv_days / 30.4
    clin$pathway = pathway

    survfit(Surv(surv_months, alive) ~ pathway, data=clin) %>%
        ggsurv() +
            xlim(0, 52) +
            theme_bw() +
            xlab("Survival (weeks)") +
            ggtitle(paste(row['subset'], row['scores'], row['adj.p'])) %>%
        print()
}

scores = io$load(INFILE) %>%
    ar$map(scores, along=1, subsets=clinical$study, discretize_quartiles)

assocs.pan = util$pancan(scores)
assocs.tissue = util$tissue(scores)

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

assocs.pan %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1) %>%
    print()

assocs.pan %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(5) %>%
    apply(1, row2survFit)

assocs.tissue %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) %>%
    print()

assocs.tissue %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(15) %>%
    apply(1, row2survFit)
