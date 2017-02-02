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
#' @param transform  transform density to normal distribution
#' @param error      return value on error (e.g. an empty list given; NULL: drop)
#' @param ...        arguments passed to wGSEA
#' @return           matrix or vector of enrichment scores (cell lines x signatures)
runGSEA = function(expr, sigs, transform.normal=FALSE, error=NA, ...) {
    if (!is.list(sigs))
        sigs = list(sigs)

    if (length(intersect(c(unlist(sigs)), rownames(expr))) == 0)
        stop("no common names between expression and signature")

    result = do.call(cbind,
        setNames(lapply(sigs, function(sig) {
            score = apply(expr, 2, function(e) wGSEA(e, sig, ...))
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
#' @param transform  transform density to normal distribution
#' @param error      return value on error (e.g. an empty list given; NULL: drop)
#' @return           matrix or vector of enrichment scores (cell lines x signatures)
runTwoTailedGSEA = function(expr, upreg, downreg, transform=TRUE, error=NA) {
    up = runGSEA(expr, upreg, transform, error)
    down = runGSEA(expr, downreg, transform, error)
    up - down
}

#' Performs weighted Gene Set Enrichment Analysis
#'
#' @param norm_express  named vector of expression values
#' @param signature     vector of genes in the target set
#' @param p             weight multiplier
#' @param display       draw enrichment plot [T/F]
#' @param returnRS      hit-miss difference for each position [T/F] {defunct}
#' @param normalize     return normalized score instead of raw
#' @param significance  compute significance by shuffling genes
#' @param trial         number of runs to shuffle genes if significance=T
#' @return              enrichtment score, or list of options specified
wGSEA = function(norm_express, signature, p=1, display=FALSE, returnRS=FALSE,
        normalize=FALSE, significance=FALSE, trial=1000) {
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
    
    if (normalize || significance) {
        EMPES = computeSimpleEMPES(norm_express, signature, trial)
        if (significance) {
#            P = ((ES>=0) * sum(EMPES>=ES) + (ES<0) * sum(EMPES<=ES)) / trial
#            return(list(NES=NES, p.value=P))
            normalizeESgamma(ES, EMPES)
        } else {
            normalizeESsimple(ES, EMPES)
        }
    } else
        ES

#    if (returnRS) {
#        POSITIONS = which(HITS == 1)
#        names(POSITIONS) = ranking[POSITIONS]
#        POSITIONS = POSITIONS[order(names(POSITIONS))]
#        return(list(ES = ES, RS = RS, POSITIONS = POSITIONS, PEAK = peak))
#    }
}

#' Uses a kernel density estimator to obtain normally distributed scores
#'
#' @param ES         a vector of enrichment scores
#' @param transform  transform the resulting distribution from Uniform to Normal [T/F]
#' @return           the resulting score
normalizeCDF = function(ES, transform.normal=TRUE) {
    CDF = sROC::kCDF(ES, xgrid=ES, adjust=1, na.rm=TRUE)
    nep = CDF$Fhat[match(ES, CDF$x)]
    if (transform.normal)
       log(nep/(1-nep))
    else
        nep
}

#' Computes a background distribution by using random gene sets
#'
#' @param exp_value_profile  named expression vector
#' @param sig     the original signature (to match lengths)
#' @param trials  how many random sets to test
#' @return        vector of enrichment scores for random sets
computeSimpleEMPES = function(exp_value_profile, sig, trials=1000) {
    sapply(1:trials, function(i) wGSEA(exp_value_profile,
        signature = sample(names(exp_value_profile), size=length(sig)),
        normalize = FALSE))
}

#' Normalize enrichment scores by background distribution
#'
#' @param ES     the raw enrichment score
#' @param empES  random empirical (as in, output of EMPES) score
#' @return       normalized enrichment score
normalizeESsimple = function(ES, empES) {
    ES[ES>0] = ES[ES>0] / mean(empES[empES>0])
    ES[ES<0] = - ES[ES<0] / mean(empES[empES<0])
    ES
}

#' Normalize enrichment score by fitting a gamma distribution
#'
#' @param ES     the raw enrichment score
#' @param empES  random empirical (as in, output of EMPES) score
#' @return       list of normalized enrichment score, p-value
normalizeESgamma = function(ES, empES) {
    if (ES >= 0){
        NULLES = fitdistr(empES[empES>=0], "gamma")
        pval = 1 - pgamma(ES, shape = NULLES$estimate[1], rate = NULLES$estimate[2])
        NES = ES / mean(empES[empES>=0])
    }
    else {
        NULLES = fitdistr(-empES[empES<0], "gamma")
        pval = 1 - pgamma(-ES, shape = NULLES$estimate[1], rate = NULLES$estimate[2])
        NES = - ES / mean(empES[empES<0])
    }
    list(NES = NES, pvalue = pval)
}

#' Takes a list of gene sets and returns the list filtered by valid IDs and number
#'
#' @param genesets  A list of vectors
#' @param valid     A vector of identifiers that can be used
#' @param min       The minimum number of genes in a list to keep the list
#' @param max       The maximum number of genes in a list to keep the list
filter_genesets = function(genesets, valid, min=5, max=500, warn=TRUE) {
	if (any(is.na(valid)))
		warning("NA found in valid set")
	if (any(valid == ""))
		warning("empty identifier found in valid set")

	genesets = lapply(genesets, function(x) intersect(x, valid))

	num_overlap = sapply(genesets, length)
	discard = num_overlap < min | num_overlap > max

	if (any(discard) && warn) {
		warning("Discarding ", sum(discard), " (of ", length(discard), ") sets: ",
                paste(names(genesets)[discard], collapse=", "))
		genesets = genesets[!discard]
	}

    if (length(genesets) == 0)
        stop("No gene sets left after filtering")

	genesets
}
