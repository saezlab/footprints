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

min_conc = gdsc$drug$conc('min', colnames(Ys))
max_conc = gdsc$drug$conc('max', colnames(Ys))

pdf(OUTFILE, width=26, height=20)

for (drug in sort(colnames(Ys))) {
    df = data.frame(tissue=tissues, drug=Ys[,drug])
    med = function(x) median(x, na.rm=TRUE)

    minc = log10(min_conc[drug])
    maxc = log10(max_conc[drug])
    rect = data.frame(xmin=-Inf, xmax=Inf, ymin=minc, ymax=maxc)

    box = ggplot(df, aes(x=reorder(tissue, drug, FUN=med), y=drug)) +
        geom_boxplot() +
        xlab("Cancer type") +
        ylab("- log IC50") +
        theme_bw() +
        geom_abline(intercept=minc, slope=0, linetype="dotted") +
        geom_abline(intercept=maxc, slope=0, linetype="dotted") +
        geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                  fill="plum2", alpha=0.1, inherit.aes=FALSE) +
        ggtitle(drug)

    print(box)
}

dev.off()
