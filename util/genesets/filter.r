#' Takes a list of gene sets and returns the list filtered by valid IDs and number
#'
#' @param genesets  A list of vectors
#' @param valid     A vector of identifiers that can be used
#' @param min       The minimum number of genes in a list to keep the list
#' @param max       The maximum number of genes in a list to keep the list
filter = function(genesets, valid, min=5, max=500, warn=TRUE) {
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
