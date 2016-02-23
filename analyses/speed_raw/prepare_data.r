library(dplyr)
library(tidyr)
b = import('base')
io = import('io')
ar = import('array')

data = io$load('../../data/expr.RData')

remove_field = function(x, field, from, to) {
    x[[field]] = NULL
    x[[to]] = x[[from]]
    x[[from]] = NULL
    x
}

subset_expr = function(id) {
    ref = strsplit(id, ":")[[1]]
    data$expr[[ref[1]]][, ref[2]] %catch% NA
}

treat = setNames(c("o", "+", "-"), c("control", "activating", "inhibiting"))

controls = data$records %>%
    lapply(function(x) remove_field(x, "perturbed", "control", "array")) %>%
    lapply(as.data.frame) %>%
    bind_rows() %>%
    mutate(treatment = "o", hours=NA, pathway=NA)

perturbed = data$records %>%
    lapply(function(x) remove_field(x, "control", "perturbed", "array")) %>%
    lapply(as.data.frame) %>%
    bind_rows() %>%
    mutate(treatment = treat[effect])

index = bind_rows(controls, perturbed) %>%
    select(-exclusion, -effect) %>%
    mutate(expr_ref = paste(id, array, sep=":"),
           new_id = paste(accession, array, sep=":")) %>%
    distinct(new_id)

expr = sapply(index$expr_ref, subset_expr, simplify=FALSE, USE.NAMES=TRUE) %>%
    ar$stack(along=2)

stopifnot(index$expr_ref == colnames(expr))
colnames(expr) = index$new_id

index = index %>%
    mutate(id = new_id) %>%
    select(-new_id, -expr_ref)

save(index, expr, file="expr.RData")
