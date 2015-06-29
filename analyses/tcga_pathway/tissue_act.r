library(ggplot2)
library(reshape2)
library(dplyr)
io = import('io')
ar = import('array')
icgc = import('data/icgc')

speed = io$load('../../model/model_norm.RData')
avail = icgc$available(rna_seq=TRUE)

STUDY = c("Breast_Invasive_Carcinoma-TCGA-US",
        # "Pancreatic_Cancer-QCMG-AU", # no normals
        # "Ovarian_Serous_Cystadenocarcinoma-TCGA-US", # no normals
          "Lung_Adenocarcinoma-TCGA-US",
          "Kidney_Renal_Clear_Cell_Carcinoma-TCGA-US",
          "Thyroid_Carcinoma-TCGA-US",
          "Neuroblastoma-TARGET-US",
          "Uterine_Corpus_Endometrioid_Carcinoma-TCGA-US",
        # "Acute_Lymphoblastic_Leukemia-TARGET-US",
          "Colon_Adenocarcinoma-TCGA-US",
          "Lung_Squamous_Cell_Carcinoma-TCGA-US",
          "Head_and_Neck_Squamous_Cell_Carcinoma-TCGA-US")

pdf("plot.pdf", paper="a4r", width=26, height=20)
on.exit(dev.off())

for (study in STUDY) {
    # get clinical + expression data
    cc = icgc$clinical() %>%
        filter(study == study,
               specimen_type %in% c("Normal - tissue adjacent to primary",
                                    "Normal - solid tissue",
                                    "Primary tumour - solid tissue"),
               donor_sex == "female",
               icgc_specimen_id %in% avail) %>%
        select(icgc_donor_id, icgc_specimen_id, specimen_type)

    # compute speed scores for all
    expr = na.omit(icgc$rna_seq(cc$icgc_specimen_id))
    ar$intersect(speed, expr, along=1)
    scores = t(expr) %*% speed %>% ar$map(along=1, scale)
    scores = scores[cc$icgc_specimen_id,]

    # see which pathways change from normal to tumour
    df = as.data.frame(scores)
    df$type = cc$specimen_type
    df = melt(df)

    ggplot(df, aes(x=variable, y=value, fill=type)) +
        geom_boxplot() +
        xlab("pathways") +
        ylab("standard deviations") +
        ggtitle(paste0(study, ": pathway activation normal vs tumour")) %>%
        print()
}
