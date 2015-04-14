#!/usr/bin/env Rscript
io = import('io')
plt = import('plots')

result = na.omit(io$load('models.RData'))
assocs = io$load('assocs.RData')

pdf("plot.pdf", paper="a4r", width=26, height=20)

### draw volcano plots
plt$drawVolcano(assocs$assocs$speed, top=40, log='y', base.size=0.2) +
    ggtitle("SPEED pathway associations that can not be explained by TCGA label, 5% FDR") + 
    xlab("Regression slope") + 
    ylab("FDR-adjusted p-value (log)")

plt$drawVolcano(assocs$assocs$opt50, top=40, log='y', base.size=0.2) +
    ggtitle("opt SPEED pathway associations that can not be explained by TCGA label, 5% FDR") + 
    xlab("Regression slope") + 
    ylab("FDR-adjusted p-value (log)")

plt$drawVolcano(assocs$assocs$gatza, top=40, log='y', base.size=0.2) +
    ggtitle("Gatza pathway associations that can not be explained by TCGA label, 5% FDR") + 
    xlab("Regression slope") + 
    ylab("FDR-adjusted p-value (log)")

### draw speed vs opt50 EN results
plt$drawModelCompare(result, "speed", "gatza") +
    ggtitle("Elastic Net cross-validation errors for <= 4 features (lower is better)")

boxplot(result[c('speed','gatza')])

### draw correlation between scores in different cell lines
library(corrplot)
corrplot(assocs$corMat$speed)
corrplot(assocs$corMat$opt50)
corrplot(assocs$corMat$gatza)

dev.off()
