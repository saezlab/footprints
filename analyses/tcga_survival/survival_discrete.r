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

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_linear.pdf"

discretize_quartiles = function(x) {
    qq = quantile(x)
    re = rep(0, length(x))
    re[x > unname(qq[4])] = 1
    re[x < unname(qq[2])] = -1
    if (sum(abs(re)) < 10)
        return(rep(NA, length(x)))

    re = as.character(re)
    re[re == "-1"] = "down"
    re[re == "0"] = "unknown"
    re[re == "1"] = "up"
    re
}

#TODO: make sure this works with util if supplying discretized scores
# and then use this instead of the code here (which is kind of crap)
row2survFit = function(row, include_normal=FALSE) {
    score = scores
    clin = util$clinical
    if ("subset" %in% names(row))
        clin = filter(util$clinical, study == row['subset'])

    ar$intersect(score, clin$barcode, along=1)
    clin$pathway = score[,sub("_.*$", "", row['scores'])]

    survfit(Surv(surv_months, alive) ~ pathway, data=clin) %>%
        ggsurv() +
            xlim(0, 52) +
            theme_bw() +
            xlab("Survival (weeks)") +
            ggtitle(paste(row['subset'], row['scores'], row['adj.p'])) %>%
        print()
}

scores = io$load(INFILE)
clinical = util$clinical

# select primary tumors only, one sample per patient
scores = scores[substr(rownames(scores), 14, 16) == "01A",]
rownames(scores) = substr(rownames(scores), 1, 12)

ar$intersect(scores, clinical$barcode, along=1)

scores = scores %>%
    ar$map(along=1, subsets=clinical$study, discretize_quartiles) %>%
    as.data.frame() %>%
    lapply(function(x) setNames(relevel(factor(x), "unknown"), rownames(.))) %>%
    do.call(data.frame, .)

assocs.pan = util$pancan(scores, clinical)
assocs.tissue = util$tissue(scores, clinical)

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

assocs.pan %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = scores) %>%
    plt$volcano(base.size=0.1) %>%
    print()

fits = assocs.pan %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(5)
if (nrow(fits) >= 1)
    apply(fits, 1, row2survFit)

assocs.tissue %>%
    plt$color$p_effect("adj.p", dir=-1) %>%
    mutate(label = paste(subset, scores, sep=":")) %>%
    plt$volcano(p=0.1) %>%
    print()

fits = assocs.tissue %>%
    arrange(adj.p) %>%
    filter(adj.p < 0.1) %>%
    head(15)
if (nrow(fits >= 1))
    apply(fits, 1, row2survFit)
