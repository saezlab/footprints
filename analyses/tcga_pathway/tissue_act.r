library(ggplot2)
library(reshape2)
library(dplyr)
io = import('io')
ar = import('array')
df = import('data_frame')
icgc = import('data/icgc')

speed_full = io$load('../../model/model_norm.RData')
avail = icgc$available(rna_seq=TRUE)

STUDY = c("Breast_Invasive_Carcinoma-TCGA-US",
        # "Pancreatic_Cancer-QCMG-AU", # no normals
        # "Ovarian_Serous_Cystadenocarcinoma-TCGA-US", # no normals
          "Lung_Adenocarcinoma-TCGA-US",
          "Kidney_Renal_Clear_Cell_Carcinoma-TCGA-US",
          "Thyroid_Carcinoma-TCGA-US",
#        "Neuroblastoma-TARGET-US", # nothing avail
          "Uterine_Corpus_Endometrioid_Carcinoma-TCGA-US",
        # "Acute_Lymphoblastic_Leukemia-TARGET-US", # blood cancer, needs different fields
          "Colon_Adenocarcinoma-TCGA-US",
          "Lung_Squamous_Cell_Carcinoma-TCGA-US",
          "Head_and_Neck_Squamous_Cell_Carcinoma-TCGA-US")

pdf("plot.pdf", paper="a4r", width=26, height=20)
on.exit(dev.off())

for (cur_study in STUDY) {
    # get clinical + expression data
    cc = icgc$clinical() %>%
        filter(study == cur_study,
               specimen_type %in% c("Normal - tissue adjacent to primary",
                                    "Normal - solid tissue",
                                    "Primary tumour - solid tissue"),
#               donor_sex == "female",
               icgc_specimen_id %in% avail) %>%
        select(icgc_donor_id, icgc_specimen_id, specimen_type, tissue)

    # compute speed scores for all
    speed = speed_full
    expr = na.omit(icgc$rna_seq(cc$icgc_specimen_id, voom=TRUE))
    ar$intersect(speed, expr, along=1)
    scores = t(expr) %*% speed %>% ar$map(along=1, scale)
    scores = scores[cc$icgc_specimen_id,]

#    # compute linear model for pathway changes
#    # because of the sample sizes everything is significant - need other method
#    design = ar$mask(cc$specimen_type) + 0
#    colnames(design) = c("tumor", "normal")
#    fit.1 = limma::lmFit(t(scores), design)
#    contrast = makeContrasts("tumor-normal", levels=design)
#    fit.2 = limma::contrasts.fit(fit.1, contrast)
#    fit.3 = limma::eBayes(fit.2)
#    pval = setNames(p.adjust(fit.3$p.value, method="fdr"), rownames(fit.3))
#    coeff = fit.3$coefficients[,1]

#    # calculate fold change and distribution overlap
#    means = ar$map(scores, along=1, mean, subsets=cc$specimen_type)
#    sds = ar$map(scores, along=1, sd, subsets=cc$specimen_type)
#    calls = data.frame(mu1 = means[1,], mu2 = means[2,],
#                       sd1 = sds[1,], sd2 = sds[2,])
#    overlaps = df$call(calls, st$overlap_normals, result_only=TRUE)

    # boxplot which pathways change from normal to tumour
    df = as.data.frame(scores)
    df$type = cc$specimen_type
    df = melt(df)

    print(ggplot(df, aes(x=variable, y=value, fill=type)) +
        geom_boxplot() +
        xlab("pathways") +
        ylab("standard deviations") +
        ggtitle(paste0(cur_study, ": pathway activation normal vs tumour")))
}
