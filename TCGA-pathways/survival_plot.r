library(limma)
library(dplyr)
b = import('base')
io = import('io')
an = import('anova')
plt = import('plots')

# load clinical data, extract survival time+status
clinical = io$load('survival.RData')
scores = io$load('expr_scores.RData')

survival = clinical[c('known_survival_time','donor_vital_status')]
rownames(survival) = clinical$icgc_sample_id

# run cox regression for associations (tissue as covariate)
assocs = an$calcAssocs(Ys = survival,
                       Xs = scores,
                       subsets = clinical$tissue,
                       type = "cox")

pdf("survival.pdf", paper="a4r", width=26, height=20)
print(plt$drawVolcano(assocs, top=40, log='y', base.size=0.2))
dev.off()
