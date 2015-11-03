library(corrplot)
library(ggplot2)
library(reshape2)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
st = import('stats')
tcga = import('data/tcga')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_linear.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "tissue_act.pdf"

# create simple data.frame with barcode, type, and study; same for scores
scores = io$load(INFILE)
index = tcga$barcode2index(rownames(scores)) %>%
    filter(grepl("Primary|Normal", Sample.Definition),
           Vial == "A") %>%
    transmute(barcode = Bio.ID,
              type = Sample.Definition,
              study = Study.Abbreviation) %>%
    mutate(type = ifelse(grepl("Normal", type), "normal", "tumor"))
scores = scores[index$barcode,]

pdf(OUTFILE, paper="a4r", width=26, height=20)
on.exit(dev.off)

# box + violin plot for each tissue
for (tissue in unique(index$study)) {
    message(tissue)
    stopifnot(rownames(scores) == index$barcode)
    score = scores[index$study == tissue,]
    type = filter(index, study==tissue)$type
    normals = type == "normal"

    if (sum(normals) < 5)
        next


    cor_normal = cor(score[normals,])
    cor_cancer = cor(score[!normals,])

    # median=0, sd=1 for the normals, same factor for tumors
    score = score %>% ar$map(along=1, function(x) {
        (x - median(x[normals])) / sd(x[normals])
    })

    df = melt(data.frame(type=type, score, check.names=FALSE))

    box = ggplot(df, aes(x=variable, y=value, fill=type)) +
        geom_boxplot() +
        xlab("pathways") +
        ylab("standard deviations") +
        ggtitle(paste0(tissue, ": pathway activation normal vs tumour"))

#    violin = ggplot(df, aes(x=variable, y=value, fill=type)) +
#        geom_violin() +
#        xlab("pathways") +
#        ylab("standard deviations") +
#        ggtitle(paste0(tissue, ": pathway activation normal vs tumour"))

    print(box)
#    print(violin)

    old.par = par(mfrow=c(1, 2))
    corrplot(cor_normal)
    corrplot(cor_cancer)
    par(old.par)
}
