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
OUTFILE = commandArgs(TRUE)[2] %or% "correlation.RData"

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

tissue2cor = function(tissue) {
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
        return(NULL)

    # calculate pathway cors in normal, tumor, and difference
    cor_normal = cor(score[normals,])
    cor_cancer = cor(score[!normals,])
    cor_diff = 0.5 * (cor_cancer - cor_normal)
    p_diff = st$cor$diff_test(score[!normals,], score[normals,])

    ar$stack(list(normal = cor_normal, cancer = cor_cancer,
                  diff = cor_diff, p.value = p_diff), along=3)
}

iter = c("pan", sort(unique(index$study)))
cor = lapply(iter, tissue2cor) %>%
    setNames(iter) %>%
    ar$stack(along=4) %>%
    melt() %>%
    tidyr::spread(Var3, value) %>%
    transmute(path1 = Var1,
              path2 = Var2,
              tissue = Var4,
              normal = normal,
              cancer = cancer,
              diff = diff,
              p.value = p.value) %>%
    filter(path1 != path2) %>%
    arrange(p.value)

cor = cor[seq(1, nrow(cor), 2),]

save(cor, file=OUTFILE)
