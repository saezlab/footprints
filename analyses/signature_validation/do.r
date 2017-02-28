library(ggplot2)
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

    re = list()

    if (!is.null(idx$expr$description))
        key = paste("PROGENy", idx$expr$description, sep=": ")
    else
        key = "PROGENy"

    re[[key]] = list(
        control = scores[idx$expr$control, path],
        perturbed = scores[idx$expr$perturbed, path]
    )
    re[[idx$activity$description]] = list(
        control = idx$activity$control,
        perturbed = idx$activity$perturbed
    )

    df = reshape2::melt(re) %>%
        transmute(pathway = paste(path, idx$paper$description, sep="\n"),
                  method = L1,
                  type = L2,
                  value = value)
}

idf = lapply(names(index), pathway_scores) %>%
    dplyr::bind_rows() %>%
    group_by(pathway, method) %>%
    mutate(value = (value - median(value[type == "control"])) /
           sd(value[type == "control"])) %>%
    ungroup() %>%
    group_by(pathway, method) %>%
    mutate(value = scale(value, center=FALSE)) %>%
    ungroup()

p = ggplot(idf, aes(x=method, y=value, fill=type)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(shape=21, size=3, alpha=0.8, position=position_dodge(0.75)) +
    facet_wrap(~pathway, scales="free") +
    theme_minimal() +
    xlab("") +
    ylab("A.U.")

if (is.null(module_name())) {
    pdf("plot.pdf")
    print(p)
}
