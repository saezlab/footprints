library(ggplot2)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

model = io$load(module_file('../../model/model_matrix.RData'))$model
expr = io$load(module_file('expr.RData'))
index = io$read_yaml(module_file('validation.yaml'))

# update index to include scores instead of expr index
pathway_scores = function(path) {
    message(path)
    idx = index[[path]]

    # calculate pathway scores
    eobj = expr[[idx$expr$accession]]
    pd = Biobase::pData(eobj)
    exp = Biobase::exprs(eobj)
    colnames(exp) = pd[colnames(exp), 'Source.Name']
    mod = model
    ar$intersect(mod, exp, along=1)
    scores = t(exp) %*% mod

    if (is.null(idx$expr$description))
        idx$expr$description = "PROGENy"
    else
        idx$expr$description = paste("PROGENy", idx$expr$description, sep=": ")

    re = list()
    re[[idx$expr$description]] = list(
        control = scores[idx$expr$control, path],
        perturbed = scores[idx$expr$perturbed, path]
    )
    re[[idx$activity$description]] = list(
        control = idx$activity$control,
        perturbed = idx$activity$perturbed
    )

    re2pval = function(e) {
        p = t.test(e$control, e$perturbed)$p.value / 2 # we know direction
        sprintf("p %.2g", p)
    }
    ps = lapply(re, function(r) re2pval(r) %catch% "")
    names(re) = paste(names(re), ps, sep="\n")

    df = reshape2::melt(re) %>%
        transmute(pathway = paste(path, idx$paper$description, sep="\n"),
                  method = L1,
                  Type = L2,
                  value = value)
}

idf = lapply(names(index), pathway_scores) %>%
    dplyr::bind_rows() %>%
    group_by(pathway, method) %>%
    mutate(value = (value - median(value[Type == "control"])) /
           sd(value[Type == "control"])) %>%
    ungroup() %>%
    group_by(pathway, method) %>%
    mutate(value = scale(value, center=FALSE)) %>%
    ungroup()

idf$method = factor(idf$method)
pidx = grepl("PROGENy", levels(idf$method))
idf$method = forcats::fct_relevel(idf$method,
    levels(idf$method)[c(which(!pidx), which(pidx))])

p = ggplot(idf, aes(x=method, y=value, fill=Type)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(shape=21, size=3, alpha=0.8, position=position_dodge(0.75)) +
    facet_wrap(~pathway, scales="free") +
    theme_minimal() +
    theme(legend.position = c(0.92, 0.1),
          legend.justification = c(1, 0),
          legend.text = element_text(size=12)) +
    xlab("") +
    ylab("A.U.")

if (is.null(module_name())) {
    pdf("plot.pdf", paper="a4r", width=26, height=20)
    print(p)
}
