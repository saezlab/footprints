library(ggplot2)
library(dplyr)
io = import('io')
ar = import('array')
gdsc = import('data/gdsc')

tissues = gdsc$tissues("SCLC")
Ys = gdsc$drug_response('IC50s')
mut = gdsc$mutated_genes()
scores = io$load('../../scores/gdsc/speed_norm.RData')
ar$intersect(Ys, mut, scores, tissues, along=1)

df = data.frame(mut = mut[,"MAP2K4"] | mut[,"MAP3K1"] |
                    mut[,"HRAS"] | mut[,"KRAS"] | mut[,"NRAS"] | mut[,"BRAF"],
                speed = scale(scores[,"MAPK"]),
                drug = Ys[,'PD-0325901'])

# takes a data.frame, and returns only matching columns with value
dd = function(df, condition, value, label, field="value") {
    sub = df[condition,]
    sub[[field]] = value
    sub$label = label
    sub
}

data = rbind(
    dd(df, !df$mut, "wt", "wtmut"),
    dd(df, df$mut, "mut", "wtmut"),
    dd(df, df$speed > 1, "active", "pathway"),
    dd(df, abs(df$speed) <= 1, "unknown", "pathway"),
    dd(df, df$speed < -1, "inactive", "pathway"),
    dd(df, df$mut & df$speed > 1, "active+mut", "both")
)

data$value = factor(data$value, levels=c("wt","mut","inactive","unknown","active","active+mut"))

# histogram: drug response for mutated/wt
#ggplot(df, aes(x=mut, y=drug)) + geom_boxplot()

pdf("SCLC-pd-mapk.pdf", width=10, height=8)
on.exit(dev.off)

ggplot(data, aes(x=value, y=drug, fill=label)) +
    geom_boxplot(outlier.size=NA) +
    geom_point(shape=21, colour="grey", position=position_jitter(width=.25)) +
    theme_bw() +
    xlab("Phenotype") +
    ylab("Drug IC50 [log uM]") +
    ggtitle("SCLC response to PD-0325901 [MAPK]")
#ggplot(df, aes(x=p53, y=drug, fill=tissue)) + geom_boxplot() + facet_grid(.~tissue)
