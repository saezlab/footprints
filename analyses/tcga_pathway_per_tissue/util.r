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
gdsc = import('data/gdsc')

#' Function to return pathway scores and index
get_tumors = function(method, tissue=NULL) {
    # create simple data.frame with barcode, type, and study; same for scores
    scores = io$load(module_file("../../scores/tcga/pathways_mapped",
                     paste0(method, ".RData")))

    index = tcga$barcode2index(rownames(scores)) %>%
        filter(grepl("Primary|Normal", Sample.Definition) & Vial == "A") %>%
        transmute(barcode = Bio.ID,
                  type = ifelse(grepl("Normal", Sample.Definition),
                                "normal", "tumor"),
                  study = Study.Abbreviation)

    c(list(scores = scores[index$barcode,]), index)
}

#' Function to return pathway scores and index
get_clines = function(method, tissue=NULL) {
    scores = io$load(module_file("../../scores/gdsc/pathways_mapped",
                     paste0(method, ".RData")))

    c(list(scores = scores), gdsc$cosmic$id2tissue(rownames(scores)))
}

#' Function to return pathway scores and index
get_merge = function(method, tissue=NULL) {
    scores = io$load(module_file("../../scores/merge/pathways_mapped",
                     paste0(method, ".RData")))

	tumors = tcga$barcode2index(rownames(scores)) %>%
		filter(grepl("Primary|Normal", Sample.Definition) & Vial == "A") %>%
        transmute(barcode = Bio.ID,
                  type = ifelse(grepl("Normal", Sample.Definition),
                                "normal", "tumor"),
                  study = Study.Abbreviation)

	clines = gdsc$tissues()[grep("^[0-9]+$", rownames(scores), value=TRUE)] %>%
		data.frame(barcode = names(.), type="cell line", study=unname(.))

	index = bind_rows(tumors, clines)
    c(scores = list(scores[index$barcode,]), index)
}

#' Compare the correlation of pathways between normal and tumor
#'
#' @param method       Method to use for comparison
#' @param tissue       Restrict analysis to cancer type (default: pan-cancer)
#' @param diff         Add a plot for the difference with -log(p-value)
#' @param min_normals  Minimum number of matched controls to compute diff
#' @return             Correlation plot object
cor_tumor_norm = function(method, tissue=NULL, diff=TRUE, min_normals=10) {
    data = get_tumors(method=method, tissue=tissue)
	score_ref = data$scores[data$type == "normal",]
	score_test = data$scores[data$type == "tumor",]

    # calculate pathway cors in normal, tumor, and difference
    cor_ref = cor(score_ref)
    cor_test = cor(score_test)
    cor_diff = 0.5 * (cor_test - cor_ref)
    p_diff = -log10(st$cor$diff_test(score_test, score_ref))
    p_diff[p_diff<2 | p_diff==Inf] = 0
    p_diff = floor(p_diff)

    old.par = par(mfrow=c(1, 3))
    corrplot(cor_ref, title="normal")
    corrplot(cor_test, title="cancer")
    corrplot(cor_diff, title="difference",
             p.mat=p_diff, sig.level=1, insig="p-value")
    par(old.par)
}

#' Compare the correlation of pathways between cell lines and tumors
#'
#' @param method  Method to use for comparison
#' @param tissue  Restrict analysis to cancer type (default: pan-cancer)
#' @param diff    Add a plot for the difference with -log(p-value)
#' @return        Correlation plot object
cor_cline_tumor = function(method, tissue=NULL, diff=TRUE) {
    tumor = get_tumors(method=method, tissue=tissue)
    clines = get_clines(method=method, tissue=tissue)
	score_ref = tumor$scores[tumor$index$type == "tumor",]
	score_test = clines$scores

    # calculate pathway cors in normal, tumor, and difference
    cor_ref = cor(score_ref)
    cor_test = cor(score_test)
    cor_diff = 0.5 * (cor_test - cor_ref)
    p_diff = -log10(st$cor$diff_test(score_test, score_ref))
    p_diff[p_diff<2 | p_diff==Inf] = 0
    p_diff = floor(p_diff)

    old.par = par(mfrow=c(1, 3))
    corrplot(cor_ref, title="tumor")
    corrplot(cor_test, title="cell line")
    corrplot(cor_diff, title="difference",
             p.mat=p_diff, sig.level=1, insig="p-value")
    par(old.par)
}

#' Box plots of relative pathway activations TCGA tumors vs normals 
#'
#' @param method  Method to use for comparison
#' @param tissue  Restrict analysis to cancer type (default: pan-cancer)
#' @return        A ggplot2 object
box_tumor_normal = function(method, tissue=NULL) {
    data = get_tumors(method=method, tissue=tissue)
    score = data$scores

    if (is.null(tissue))
        assocs = filter(st$lm(scores ~ study + type, data=data), term=="typetumor")
    else
        assocs = st$lm(scores ~ type, data=data)

    assocs = assocs %>%
        transmute(variable = score,
                  value = 5, #floor(max(score, na.rm=TRUE)),
                  label = ifelse(p.value>0.05, "", format(p.value, digits=2)))

    # median=0, sd=1 for the normals, same factor for tumors
    score = score %>% ar$map(along=1, function(x) {
        (x - median(x[normals])) / sd(x[normals])
    })

    df = melt(data.frame(type=type, score, check.names=FALSE))

    ggplot(df, aes(x=variable, y=value)) +
        geom_boxplot(aes(fill=type), outlier.shape=NA) +
        xlab("pathways") +
        ylab("standard deviations") +
        theme_bw() +
        ggtitle(paste0(tissue, ": pathway activation ", sum(normals),
                       " normals vs ", sum(!normals), " tumours")) +
        geom_text(data=assocs, aes(label=label), size=3)
}

#' Draw box plots using merged tumor and cell line data
box_merge = function(method, tissue=NULL) {
	data = get_merge(method=method, tissue=tissue)
    score = data$scores
    type = data$type

    # median=0, sd=1 for the normals, same factor for tumors
    score = score %>% ar$map(along=1, function(x) {
        (x - median(x[normals])) / sd(x[normals])
    })

    df = data.frame(type=factor(type, levels=c("normal", "tumor", "cell line")),
                    score, check.names=FALSE) %>% melt()

    box = ggplot(df, aes(x=variable, y=value, fill=type)) +
        geom_boxplot() +
        xlab("pathways") +
        ylab("standard deviations") +
        ggtitle(paste0(tissue, ": pathway activation normal vs tumour")) +
        theme_bw()
}

#' Heatmap of pathway scores for TCGA normals
matrix_normals = function(method) {
	data = get_tumor(method=method, tissue=NULL)

	assocs = st$lm(scores ~ study, data=data) %>%
		mutate(tissue = sub("^study", "", term))

	mat = ar$construct(estimate ~ tissue + scores, data=assocs)
    pheatmap::pheatmap(mat)
}

#' Heatmap of pathway scores of TCGA tumors relative to normals
matrix_tumors = function() {
	data = get_tumor(method=method, tissue=NULL)

	assocs = st$lm(scores ~ type * study, data=data) %>%
        filter(grepl("typetumor:study", term)) %>%
        mutate(term = sub("typetumor:study", "", term))

	mat = ar$construct(statistic ~ tissue + scores, data=assocs)
    pheatmap::pheatmap(mat)
}
