library(dplyr)
library(magrittr)
library(bhtsneR)
library(impute)
library(limma)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% "../../data/expr.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "tsne_fc_sig.RData"

# calculate scores from expr and speed vectors
speed = io$load(INFILE)
index = speed$records
expr = speed$expr

# scaling: assume mean/sd across scores per sample is constant
# this protects against missing genes, etc in platform
expr2scores = function(index, expr) {
    mat = expr[,c(index$control, index$perturbed)]
    design = c(rep("control", length(index$control)),
               rep("perturbed", length(index$perturbed)))
    design_matrix = model.matrix(~ 0 + as.factor(design))
    colnames(design_matrix) = c("control", "perturbed")
    fit = limma::lmFit(mat, design_matrix)
    contrast = limma::makeContrasts("perturbed-control", levels=design_matrix)
    fit2 = limma::contrasts.fit(fit, contrast)
    fit3 = limma::ebayes(fit2)

    result = data.frame(hgnc = rownames(fit2),
                        fc = fit2$Amean,
                        p.value = fit3$p.value[,1]) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        filter(adj.p < 0.1) %$%
        setNames(fc, hgnc)

    if (index$effect == "activating")
        result
    else
        -result
}

scores = mapply(expr2scores, index=index, expr=expr, SIMPLIFY=FALSE) %>%
    ar$stack(along=1, fill=0)

index = index[rownames(scores)]
index = lapply(index, function(x) x[!names(x) %in% c("control","perturbed","exclusion")]) %>%
	do.call(bind_rows, .)

dim2 = tsne(scores)
index$x = dim2[,1]
index$y = dim2[,2]

save(index, file=OUTFILE)
