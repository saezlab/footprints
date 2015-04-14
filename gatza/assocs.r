io = import('io')
ar = import('array')
an = import('anova')
sg = import('sanger_robject')

speed = io$load('../GDSC-Scores/speed.reportv1.RData')$paper
opt50 = io$load('../GDSC-Scores/speed.reportv1.RData')$opt50
gatza = io$load('../GDSC-Scores/gatza.RData')[[1]]
Ys = sg$getDrugResponseForCellLines('IC50s') # or AUC
tissues = sg$getTissues(minN=5)
ar$intersect(speed, opt50, gatza, tissues, Ys, along=1)

assocs = list(speed = an$calcAssocs(Ys, speed, covariate=tissues, p.adjust="fdr"),
              opt50 = an$calcAssocs(Ys, opt50, covariate=tissues, p.adjust="fdr"),
              gatza = an$calcAssocs(Ys, gatza, covariate=tissues, p.adjust="fdr"))

corMat = list(speed = cor(speed),
              opt50 = cor(opt50),
              gatza = cor(gatza))

save(assocs, corMat, file="assocs.RData")
