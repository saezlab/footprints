library(dplyr)
io = import('io')
ar = import('array')
gsva = import('../../scores/speed/gsva')$gsva

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
    top100 = as.data.frame(mod$p.value) %>%
        mutate(gene = rownames(.)) %>%
        arrange(`typeptb-typectl`) %>%
        head(100) %$%
        gene
}

if (is.null(module_name())) {
    # expr$expr    : list[experiments] of expression matrices
    # expr$records : index as list
    expr = io$load('../../data/expr.RData')
    signatures = mapply(exp2sig, expr$expr, expr$records, SIMPLIFY=FALSE)

    # each element = one experiment, all signatures -> stack along 1
    scores = clustermq::Q(gsva, index=expr$records, expr=expr$expr,
            const = list(sigs=signatures), memory = 1024, n_jobs = 10) %>%
        setNames(names(expr$records)) %>%
        ar$stack(along=1)

    stopifnot(all(rownames(scores) == colnames(scores)))
    diag(scores) = NA

    index = lapply(expr$records, function(x) x[! names(x) %in% c('control', 'perturbed')])
    index = bind_rows(index) %>%
        select(-exclusion)

    # "pathways" are in cols
    save(index, scores, file="sigs_gsva.RData")
}
