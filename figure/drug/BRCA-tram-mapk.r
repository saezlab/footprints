library(ggplot2)
library(dplyr)
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')

tissues = gdsc$tissues("BRCA")
Ys = gdsc$drug_response('IC50s')
mut = gdsc$mutated_genes()
scores = io$load('../scores/gdsc/speed_norm.RData')
ar$intersect(Ys, mut, scores, tissues, along=1)

df = data.frame(mut = mut[,"MAP2K4"] | mut[,"MAP3K1"] |
                    mut[,"HRAS"] | mut[,"KRAS"] | mut[,"NRAS"] | mut[,"BRAF"],
                speed = scale(scores[,"MAPK"]),
                drug = Ys[,'Trametinib'])

# takes a data.frame, and returns only matching columns with value
dd = function(df, condition, value, field="value") {
    sub = df[condition,]
    sub[[field]] = value
    sub
}

data = rbind(
    dd(df, !df$mut, "wt"),
    dd(df, df$mut, "mut"),
    dd(df, df$speed > 1, "active"),
    dd(df, abs(df$speed) <= 1, "unknown"),
    dd(df, df$speed < -1, "inactive"),
    dd(df, df$mut & df$speed > 1, "active+mut")
)

data$value = factor(data$value, levels=c("wt","mut","inactive","unknown","active","active+mut"))

# histogram: drug response for mutated/wt
#ggplot(df, aes(x=mut, y=drug)) + geom_boxplot()

ggplot(data, aes(x=value, y=drug)) +
    geom_boxplot(outlier.size=NA) +
    geom_point(shape=21, position=position_jitter(width=.25)) +
    xlab("Phenotype") +
    ylab("Drug IC50 [log uM]") +
    ggtitle("BRCA response to Trametinib [MAPK]")
#ggplot(df, aes(x=p53, y=drug, fill=tissue)) + geom_boxplot() + facet_grid(.~tissue)
