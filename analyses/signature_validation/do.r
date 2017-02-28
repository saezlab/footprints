library(dplyr)
io = import('io')
ar = import('array')

model = io$load('../../model/model_matrix.RData')$model
expr = io$load('expr.RData')
index = io$read_yaml('validation.yaml')

# update index to include scores instead of expr index
pathway_scores = function(path) {
    idx = index[[path]]

    # calculate pathway scores
    eobj = expr[[idx$expr$accession]]
    pd = Biobase::pData(eobj)
    exp = Biobase::exprs(eobj)
    colnames(exp) = pd[colnames(exp), 'Source.Name']
    mod = model
    ar$intersect(mod, exp, along=1)
    scores = t(exp) %*% mod

    re = list(PROGENy = list(
        control = scores[idx$expr$control, path],
        perturbed = scores[idx$expr$perturbed, path]
    ))
    re[[idx$paper$description]] = list(
        control = idx$activity$control,
        perturbed = idx$activity$perturbed
    )

    df = reshape2::melt(re) %>%
        transmute(pathway = path,
                  method = L1,
                  type = L2,
                  value = value)
}

idf = lapply(names(index), pathway_scores) %>%
    dplyr::bind_rows() %>%
    group_by(pathway, method) %>% #TODO: scale by control
    mutate(value = (value - median(value[type == "control"])) /
           sd(value[type == "control"])) %>%
    ungroup()

# plot
