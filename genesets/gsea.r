###
### utility functions for GSEA analysis
###
library(MASS)
library(sROC)
library(parallel)

#' Calculates raw GSEA scores for a list of signatures
#'
#' @param expr       expression matrix (genes x samples)
#' @param sigs       signature gene symbol vector or list thereof
#' @param tcenter    tissue-center expression (default: false) [not implemented]
#' @param transform  transform density to normal distribution
#' @param error      return value on error (e.g. an empty list given; NULL: drop)
#' @return            matrix or vector of enrichment scores (cell lines x signatures)
runGSEA = function(expr, sigs, tcenter=F, transform.normal=F, error=NA) {
    if (!is.list(sigs))
        sigs = list(sigs)

    if (length(intersect(c(unlist(sigs)), rownames(expr))) == 0)
        stop("no common names between expression and signature")

    result = do.call(cbind,
        setNames(lapply(sigs, function(sig) {
            score = apply(expr, 2, function(e) wGSEA(e, sig))
            if (all(is.na(score)))
                error
            else {
                if (transform.normal)
                    normalizeCDF(score, TRUE)
                else
                    score
            }
        }), names(sigs))
    )

    if (!identical(unname(drop(result)), error))
        rownames(result) = colnames(expr)

    drop(result)
}

#' Calculates two-tailed GSEA scores for a list of signatures
#'
#' @param expr       expression matrix (genes x cell lines)
#' @param upreg      signature gene symbol vector of upregulated genes or list thereof
#' @param downreg    signature gene symbol vector of downregulated genes or list thereof
#' @param tcenter    tissue-center expression (default: false) [not implemented]
#' @param transform  transform density to normal distribution
#' @param error      return value on error (e.g. an empty list given; NULL: drop)
#' @return           matrix or vector of enrichment scores (cell lines x signatures)
runTwoTailedGSEA = function(expr, upreg, downreg, tcenter=F, transform=T, error=NA) {
    up = runGSEA(expr, upreg, tcenter, transform, error)
    down = runGSEA(expr, downreg, tcenter, transform, error)
    up - down
}

#' Performs weighted Gene Set Enrichment Analysis
#'
#' @param norm_express  named vector of expression values
#' @param signature     vector of genes in the target set
#' @param p             weight multiplier
#' @param display       draw enrichment plot [T/F]
#' @param returnRS      hit-miss difference for each position [T/F]
#' @param significance  compute significance by shuffling genes
#' @param trial         number of runs to shuffle genes if significance=T
#' @return              enrichtment score, or list of options specified
wGSEA = function(norm_express, signature, p=1, display=F, returnRS=F, significance=F, trial=1000) {
    if (is.null(names(norm_express)) || !is.vector(norm_express) || is.list(norm_express))
        stop("norm_express must be a named vector")

    if (!is.character(signature) || !is.vector(signature) || is.list(signature))
        stop("signature must be a character vector")

    norm_express = sort(norm_express, decreasing=T, method="shell")
    ranking = names(norm_express)
    signature = unique(intersect(signature, ranking))
    if (length(signature) == 0)
        return(NA)

    HITS = is.element(ranking, signature) + 0
    R = norm_express * HITS
    N = length(ranking)
    
    Phit = cumsum(abs(R)^p) / sum(abs(R)^p)
    Pmiss = cumsum(1-HITS) / (N-length(signature))
    RS = Phit - Pmiss
    Pdelta = abs(RS)

    peak = which(Pdelta > max(Pdelta)-.Machine$double.eps)[1]
    ES = unname(Phit[peak] - Pmiss[peak])

    if (display) {
        c = ifelse(ES>=0, "red", "green")
        plot(0:N, c(0, Phit-Pmiss), col = c, type = "l", xlim = c(0, N),
             ylim = c(-(abs(ES)+0.5*(abs(ES))), abs(ES)+0.5*(abs(ES))),
             xaxs = "i", bty = "l", axes = FALSE,
             xlab = "Gene Rank Position", ylab = "Running Sum")
        par(new = TRUE)
        plot(0:N, rep(0, N+1), col = 'gray', type = "l", new = FALSE, xlab = "", ylab = "",
             ylim = c(-(abs(ES)+0.5*(abs(ES))), abs(ES)+0.5*(abs(ES))))
        axis(side = 2)
    }    
    
    if (significance) {
        EMPES = computeSimpleEMPES(norm_express, signature, trial)
        P = ((ES>=0) * sum(EMPES>=ES) + (ES<0) * sum(EMPES<=ES)) / trial
    }
    if (returnRS) {
        POSITIONS = which(HITS == 1)
        names(POSITIONS) = ranking[POSITIONS]
        POSITIONS = POSITIONS[order(names(POSITIONS))]
        return(list(ES = ES, RS = RS, POSITIONS = POSITIONS, PEAK = peak))
    }
    if (significance) 
        return(list(ES = ES, P = P))
    else
        return(ES)
}

#' Uses a kernel density estimator to obtain normally distributed scores
#'
#' @param ES         a vector of enrichment scores
#' @param transform  transform the resulting distribution from Uniform to Normal [T/F]
#' @return           the resulting score
normalizeCDF = function(ES, transform.normal=T) {
    CDF = kCDF(ES, xgrid=ES, adjust=1)
    nep = CDF$Fhat[match(ES, CDF$x)]
    if (transform.normal)
       log(nep/(1-nep))
    else
        nep
}
