b = import('base')
io = import('io')
ar = import('array')

# similar to ../perturbation_recall/sigs_gsva.r
exp2sig = function(expr, index) {
    message(index$id)
    library(dplyr)
    library(magrittr)

    # use only available arrays
    ctlid = intersect(index$control, colnames(expr))
    ptbid = intersect(index$perturbed, colnames(expr))

    # subset control and perturbed matrices
    emat = expr[,c(ctlid, ptbid)]
    ctl = emat[,ctlid]
    ptb = emat[,ptbid,drop=FALSE]

    # get top100 DE genes
    type = c(rep("ctl", ncol(ctl)), rep("ptb", ncol(ptb)))
    design = model.matrix(~ 0 + type)
    mod = limma::lmFit(emat, design)
    contrast = limma::makeContrasts("typeptb-typectl", levels=design)
    mod = limma::contrasts.fit(mod, contrast)
    mod = limma::eBayes(mod)
    top100 = as.data.frame(-mod$t) %>%
        mutate(gene = rownames(.)) %>%
        arrange(`typeptb-typectl`) %>%
        head(100) %$%
        gene

	mod$coefficients[top100,]
}

EXPR = commandArgs(TRUE)[1] %or% "../../data/expr.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "top100_de.RData"

expr = io$load(EXPR)
sigs = mapply(exp2sig, expr$expr, expr$records, SIMPLIFY=FALSE)

mat = ar$stack(sigs, along=2, fill=0)

save(mat, file=OUTFILE)
