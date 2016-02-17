library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

# for each contrast, compute z-scores
expr2control = function(rec, emat) {
    control = emat[,rec$control, drop=FALSE]
    colnames(control) = NULL
    control
}

control2model = function(rec, emat) {
    # build speed models
    mean_control= apply(control, 1, mean)
    sd_control = apply(control, 1, sd)
    mean_perturbed = apply(perturbed, 1, mean)
    logFC = mean_perturbed - mean_control
    
    loess(sd_control ~ mean_control)
}

expr2zscores = function(data) {
    records = data$records
    expr = data$expr

    controls = mapply(expr2control, r=records, emat=expr, SIMPLIFY=FALSE) %>%
        ar$stack(along=2)
    models = apply(controls, 2, control2model

    lapply(result, function(x) x$zscores) %>%
        ar$stack(along=2)
}

if (is.null(module_name())) {
    INFILE = commandArgs(TRUE)[1] %or% "expr.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "scores.RData"

    data = io$load(INFILE)
    zscores = expr2zscores(data)

    # separate index file w/ metadata derived from yaml [preferred?]
    save(zscores, file=OUTFILE)
}
