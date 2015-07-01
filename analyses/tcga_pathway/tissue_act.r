library(ggplot2)
library(reshape2)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
st = import('stats')
icgc = import('data/icgc')

process_study = function(cur_study) {
    # get clinical + expression data
    cc = icgc$clinical() %>%
        filter(study == cur_study,
               specimen_type %in% c("Normal - tissue adjacent to primary",
#                                    "Normal - solid tissue",
                                    "Primary tumour - solid tissue"),
#               donor_sex == "female",
               icgc_specimen_id %in% avail) %>%
        select(icgc_donor_id, icgc_specimen_id, specimen_type, tissue)

    # compute speed scores for all
    speed = speed_full
    expr = na.omit(icgc$rna_seq(cc$icgc_specimen_id, voom=TRUE))
    ar$intersect(speed, expr, along=1)
    scores = t(expr) %*% speed
    scores = scores[cc$icgc_specimen_id,]

    # median=0, sd=1 for the normals, same factor for tumors
    scores = scores %>% ar$map(along=1, function(x) {
        normal = x[grepl("^Normal", cc$specimen_type)]
        (x - median(normal)) / sd(normal)
    })

    # boxplot which pathways change from normal to tumour
    df = as.data.frame(scores)
    df$type = cc$specimen_type
    df = melt(df)

    list(
        boxplot = ggplot(df, aes(x=variable, y=value, fill=type)) +
            geom_boxplot() +
            xlab("pathways") +
            ylab("standard deviations") +
            ggtitle(paste0(cur_study, ": pathway activation normal vs tumour")),

        violinplot = ggplot(df, aes(x=variable, y=value, fill=type)) +
            geom_violin() +
            xlab("pathways") +
            ylab("standard deviations") +
            ggtitle(paste0(cur_study, ": pathway activation normal vs tumour")),

        tumor_medians = ar$map(scores[cc$specimen_type == "Primary tumour - solid tissue",],
                along=1, median)
    )
}

if (is.null(module_name())) {
    speed_full = io$load('../../model/model_norm.RData')
    avail = icgc$available(rna_seq=TRUE)

    STUDY = c(BRCA = "Breast_Invasive_Carcinoma-TCGA-US",
            # "Pancreatic_Cancer-QCMG-AU", # no normals
            # "Ovarian_Serous_Cystadenocarcinoma-TCGA-US", # no normals
              LUAD = "Lung_Adenocarcinoma-TCGA-US",
              KRCC = "Kidney_Renal_Clear_Cell_Carcinoma-TCGA-US",
              THCA = "Thyroid_Carcinoma-TCGA-US",
#        "Neuroblastoma-TARGET-US", # nothing avail
              UCEC = "Uterine_Corpus_Endometrioid_Carcinoma-TCGA-US",
            # "Acute_Lymphoblastic_Leukemia-TARGET-US", # blood cancer, needs different fields
              COAD = "Colon_Adenocarcinoma-TCGA-US",
              LUSC = "Lung_Squamous_Cell_Carcinoma-TCGA-US",
              HNSC = "Head_and_Neck_Squamous_Cell_Carcinoma-TCGA-US")

    result = sapply(STUDY, process_study, simplify=FALSE, USE.NAMES=TRUE)

    pdf("tissue_act.pdf", paper="a4r", width=26, height=20)
    on.exit(dev.off())

    medians = ar$stack(lapply(result, r -> r$tumor_medians), along=2)
    pheatmap::pheatmap(medians, scale="row")
    #TODO: could be interesting to have number of clinical drugs targeting each pathway

    lapply(result, r -> print(r$boxplot))
    lapply(result, r -> print(r$violinplot))
}
