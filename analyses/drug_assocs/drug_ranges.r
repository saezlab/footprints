library(dplyr)
library(ggplot2)
library(reshape2)
b = import('base')
ar = import('array')
gdsc = import('data/gdsc')

OUTFILE = commandArgs(TRUE)[1] %or% 'drug_ranges.pdf'

tissues = gdsc$tissues()
Ys = gdsc$drug_response('IC50s')
ar$intersect(tissues, Ys, along=1)

pdf(OUTFILE, width=26, height=20)

for (drug in sort(colnames(Ys))) {
    df = data.frame(tissue=tissues, drug=Ys[,drug])
    med = function(x) median(x, na.rm=TRUE)

    box = ggplot(df, aes(x=reorder(tissue, drug, FUN=med), y=drug)) +
        geom_boxplot() +
        xlab("Cancer type") +
        ylab("- log IC50") +
        theme_bw() +
        ggtitle(drug)

    print(box)
}

dev.off()
