.io = import('io')

#' Calculates the SPIA scores given sample- and control experiments
#'
#' @param samples     Matrix of investigated experiments [genes x samples] or
#'                    character vector specifying columns in `data`
#' @param controls    Matrix of control experiments [genes x samples] or
#'                    character vector specifying columns in `data`
#' @param data        Matrix of samples and controls if the former are character vectors
#' @param organism    KEGG code for organism (default: human, 'hsa')
#' @param pathids     Character vector of KEGG pathway IDs to use (default: all)
#' @param bootstraps  Number of bootstraps to be performed (default: 2000)
#' @param p_filter    Maximum p-value (FDR) for which differential expression is considered
#' @return            A data.frame listing pathways and scores
spia = function(samples, controls, data=NULL, organism="hsa",
        pathids=NULL, bootstraps=2000, p_filter=1) {
    b = import('base')
    dplyr = import_package('dplyr')
    limma = import_package('limma')
    spia_pkg = import_package('SPIA')

    if (is.null(data)) {
        samples = as.matrix(samples)
        controls = as.matrix(controls)
        data = cbind(samples, controls)
        types = c(rep("sample", ncol(samples)), rep("control", ncol(controls)))
    } else {
        data = data[,c(samples, controls), drop=FALSE]
        types = c(rep("sample", length(samples)), rep("control", length(controls)))
    }

    # compute differential expression between tumor and normal
    design = model.matrix(~ 0 + as.factor(types))
    colnames(design) = sub("as\\.factor\\(types\\)", "", colnames(design))
    contrast = limma$makeContrasts("sample-control", levels=design)
    fit = limma$lmFit(data, design) %>%
        limma$contrasts.fit(contrast) %>%
        limma$eBayes()
    result = dplyr$data_frame(gene = rownames(fit),
        p.value = fit$p.value[,1],
        fold_change = fit$coefficients[,1]) %>%
        dplyr$mutate(adj.p = p.adjust(p.value, method="fdr"))

    # calculate SPIA scores (relevant fields: Name [KEGG name], ID [KEGG ID], tA [score])
    re = spia_pkg$spia(
        de = setNames(result$fold_change, result$gene)[result$adj.p < p_filter],
        all = result$gene,
        organism = organism,
        pathids = pathids,
        nB = bootstraps,
        plots = FALSE,
        verbose = FALSE
    ) %catch% data.frame()

#        data.frame(id=re$ID, name=re$Name, score=re$tA)
    setNames(re$tA, re$ID)
}


#' Calculates the SPIA scores given sample- and control experiments for each sample
#'
#' @param samples     Matrix of investigated experiments [genes x samples] or
#'                    character vector specifying columns in `data`
#' @param controls    Matrix of control experiments [genes x samples] or
#'                    character vector specifying columns in `data`
#' @param data        Matrix of samples and controls if the former are character vectors
#' @param organism    KEGG code for organism (default: human, 'hsa')
#' @param pathids     Character vector of KEGG pathway IDs to use (default: all)
#' @param bootstraps  Number of bootstraps to be performed (default: 2000)
#' @param p_filter    Maximum p-value (FDR) for which differential expression is considered
#' @return            A data.frame listing pathways and scores
spia_per_sample = function(samples, controls, data=NULL, organism="hsa",
        pathids=NULL, bootstraps=2000, p_filter=1) {
    b = import('base')
    ar = import('array')

# better:
#    samples = ar$split(samples, along=-1)
# this should work for both character vectors and 
# just needs to support "last dimension" (-1) for splitting
    if (is.null(data)) {
        stopifnot(is.matrix(samples) && is.matrix(controls))
        data = cbind(samples, controls)
        samples = colnames(samples)
        controls = colnames(controls)
    }

    sapply(samples, spia, controls=controls, data=data,
           simplify=FALSE, USE.NAMES=TRUE) %>%
        ar$stack(along=1)
}


#' Maps row names from HGNC symbols to Entrez Gene IDs
#'
#' @param  expr  The expression matrix (genes x samples)
#' @return       An expression matrix with Entrez Gene IDs
map_entrez = function(expr) {
    # map gene expression from HGNC to Entrez IDs
    lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
        biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"),
        filter="hgnc_symbol", values=rownames(expr), mart=.)
    rownames(expr) = lookup$entrezgene[match(rownames(expr), lookup$hgnc_symbol)]
    expr = limma::avereps(expr[!is.na(rownames(expr)),])
}


#' Returns all pathway IDs
#'
#' @param organism  The KEGG organism ID (default: "hsa" for human)
#' @return          A character vector with pathway IDs for a given organism
pathids = function(organism="hsa") {
    fname = paste(system.file("extdata", package = "SPIA"),
                  paste0("/", organism, "SPIA"), ".RData", sep = "")

    paths = .io$load(fname)
    setNames(names(paths), sapply(paths, function(x) x$title))
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
