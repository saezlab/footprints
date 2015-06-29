library(ggplot2)
library(dplyr)
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')

tissues = gdsc$tissues(minN=10)
Ys = gdsc$drug_response('IC50s')
mut = gdsc$mutated_genes()
scores = io$load('../scores/gdsc/speed_norm.RData')
ar$intersect(Ys, mut, scores, tissues, along=1)

df = data.frame(mut = mut[,"BRAF"],
                speed = scale(scores[,"MAPK"]),
                drug = Ys[,'Dabrafenib'],
                tissue = tissues)

# histogram: drug response for mutated/wt
#ggplot(df, aes(x=mut, y=drug)) + geom_boxplot()

# histogram: split wt in high/low speed score
df = df
df$p53 = "wt"
df$p53[df$mut] = "mutated" 

df2 = df
df2$p53[df$mut & df$speed > 0] = "mut_active"
df2$p53[df$mut & df$speed < 0] = "mut_inactive"
df2 = df2[! df2$p53 %in% c("wt","mutated"),]

df = rbind(df, df2)

ggplot(df, aes(x=p53, y=drug)) + geom_boxplot()
#ggplot(df, aes(x=p53, y=drug, fill=tissue)) + geom_boxplot() + facet_grid(.~tissue)
