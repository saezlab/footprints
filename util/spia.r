#' Calculates the SPIA scores given matrices of sample- and control experiments
#'
#' @param samples     Matrix of investigated experiments [genes x samples]
#' @param controls    Matrix of control experiments [genes x samples]
#' @param per_sample  Return scores for each sample instead of one combined score
#' @param organism    KEGG code for organism (default: human, 'hsa')
#' @param pathids     Character vector of KEGG pathway IDs to use (default: all)
#' @param bootstraps  Number of bootstraps to be performed (default: 2000)
#' @param p_filter    Maximum p-value (FDR) for which differential expression is considered
#' @param verbose     Echo sample(s) currently working on
#' @return            A data.frame listing pathways and scores
spia = function(samples, controls, per_sample=FALSE, organism="hsa", pathids=NULL, 
                bootstraps=2000, p_filter=1, verbose=FALSE) {
    sample2score = function(sample) {
        if (verbose)
            message(sample)
        data2 = cbind(data[,sample,drop=FALSE], data[,colnames(controls)])
        types = c("sample", rep("control", ncol(controls)))

        # compute differential expression between tumor and normal
        design = model.matrix(~ 0 + as.factor(types))
        colnames(design) = sub("as\\.factor\\(types\\)", "", colnames(design))
        contrast = limma$makeContrasts("sample-control", levels=design)
        fit = limma$lmFit(data2, design) %>%
            limma$contrasts.fit(contrast) %>%
            limma$eBayes()
        result = dplyr$data_frame(gene = rownames(fit),
            p.value = fit$p.value[,1],
            fold_change = fit$coefficients[,1]) %>%
            dplyr$mutate(adj.p = p.adjust(p.value, method="fdr"))

        # calculate SPIA scores (relevant fields: Name [KEGG name], ID [KEGG ID], tA [score])
        re = spia$spia(
            de = setNames(result$fold_change, result$gene)[result$adj.p < p_filter],
            all = result$gene,
            organism = organism,
            pathids = pathids,
            nB = bootstraps,
            plots = FALSE,
            verbose = FALSE
        ) %catch% data.frame()

        data.frame(id=re$ID, name=re$Name, score=re$tA)
    }

    b = import('base')
    dplyr = import_package('dplyr')
    spia = import_package('SPIA')
    limma = import_package('limma')
    data = cbind(samples, controls)

    if (per_sample)
        sapply(colnames(samples), sample2score, simplify=FALSE, USE.NAMES=TRUE) %>%
            dplyr$bind_rows()
    else
        sample2score(colnames(samples))
}

#' Mapping from SPEED to KEGG pathways
speed2kegg = c(
    MAPK = "04010", # MAPK signaling pathway
    EGFR = "04012", # ErbB signaling pathway
#    PPAR = "03320", # PPAR signaling pathway
    PI3K = "04150", # mTOR signaling pathway
    Wnt = "04310", # Wnt signaling pathway
#    Notch = "04330", # Notch signaling pathway
#    Insulin = "04910", # Insulin signaling pathway
    NFkB = "04064", # NF-kappa B signaling pathway
    VEGF = "04370", # VEGF signaling pathway
    `JAK-STAT` = "04630", # Jak-STAT signaling pathway
    TGFb = "04350", # TGF-beta signaling pathway
    Trail = "04210" # Apoptosis
)

#' Mapping from KEGG to SPEED pathways
kegg2speed = setNames(names(speed2kegg), unname(speed2kegg))
