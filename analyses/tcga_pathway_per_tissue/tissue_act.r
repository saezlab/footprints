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

INFILE = commandArgs(TRUE)[1] %or% "../../scores/tcga/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "box+cor/speed_matrix.pdf"

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
for (tissue in c("pan", sort(unique(index$study)))) {
    message(tissue)
    stopifnot(rownames(scores) == index$barcode)

    if (tissue == "pan") {
        score = scores
        type = index$type
    } else {
        score = scores[index$study == tissue,]
        type = filter(index, study==tissue)$type
    }
    normals = type == "normal"

    if (sum(normals) < 5)
        next

    # reformat stats summary df
    if (tissue == "pan") {
        assocs = st$lm(score ~ study + type, data=index) %>%
            filter(term == "typetumor")
    } else {
        assocs = st$lm(score ~ type)
    }

    assocs = assocs %>%
        transmute(variable = score,
                  value = 5, #floor(max(score, na.rm=TRUE)),
                  label = ifelse(p.value>0.05, "", format(p.value, digits=2)))

    # calculate pathway cors in normal, tumor, and difference
    cor_normal = cor(score[normals,])
    cor_cancer = cor(score[!normals,])
    cor_diff = 0.5 * (cor_cancer - cor_normal)
    p_diff = -log10(st$cor$diff_test(score[!normals,], score[normals,]))
    p_diff[p_diff<2 | p_diff==Inf] = 0
    p_diff = floor(p_diff)

    # median=0, sd=1 for the normals, same factor for tumors
    score = score %>% ar$map(along=1, function(x) {
        (x - median(x[normals])) / sd(x[normals])
    })

    df = melt(data.frame(type=type, score, check.names=FALSE))

    box = ggplot(df, aes(x=variable, y=value)) +
        geom_boxplot(aes(fill=type), outlier.shape=NA) +
        xlab("pathways") +
        ylab("standard deviations") +
        theme_bw() +
        ggtitle(paste0(tissue, ": pathway activation ", sum(normals),
                       " normals vs ", sum(!normals), " tumours")) +
        geom_text(data=assocs, aes(label=label), size=3)

    print(box)

    old.par = par(mfrow=c(1, 3))
    corrplot(cor_normal, title="normal")
    corrplot(cor_cancer, title="cancer")
    corrplot(cor_diff, title="difference",
             p.mat=p_diff, sig.level=1, insig="p-value") %catch% NULL
    par(old.par)
}
