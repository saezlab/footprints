library(magrittr)
library(dplyr)
io = import('io')
ar = import('array')
df = import('data_frame')
gdsc = import('data/gdsc')
vp = import('../analyses/drug_assocs/plot')

tissues = gdsc$tissues()
Ys = gdsc$drug_response('IC50s')
scores = io$load("../scores/gdsc/pathways_mapped/speed_matrix.RData")
mut = gdsc$mutated_genes(intogen=TRUE)
ar$intersect(tissues, Ys, scores, mut, along=1)

#' Volcano plots
volcano = function(fid) {
    fp = io$file_path('assocs_mapped', fid, ext=".RData")
    assocs = vp$load_fun(fp)$assocs.pan
    vp$plot_pancan(assocs, text.size=2)
}

#' Returns names of a logical vector where the elements are TRUE
nst = function(x) names(x[x])

#' Generates a stratification structure using pathway and mutations
#'
#' @param pathway  Our pathway descriptor, e.g. "MAPK"
#' @param genes    A character vector of which genes should be
#'                 used as markers to compare to
#' @return         A list of subsets each with a character of COSMIC IDs
stratify = function(pathway="MAPK", genes=c("BRAF","KRAS","NRAS")) {
    path = scores[names(tissues), pathway]
    has_mut = apply(mut[, genes, drop=FALSE], 1, any)
    path_active = path > quantile(path)[4] # top 25%
    path_inactive = path < quantile(path)[2] # bottom 25%
    path_null = path > quantile(path)[2] & path < quantile(path)[4] # rest

    re = list(
        # stratification by mutation
        "GENE_all" = nst(!is.na(path)),
        "GENE_wt" = nst(!has_mut),
        "GENE_mut" = nst(has_mut),

        # stratification by SPEED score
        "PATH+" = nst(path_active),
        "PATH-" = nst(path_inactive),
        "PATH_null" = nst(path_null),

        # stratification by SPEED score, wild-type subset
        "PATH+_wt" = nst(!has_mut & path_active),
        "PATH-_wt" = nst(!has_mut & path_inactive),

        # stratification by SPEED score, mutated subset
        "PATH+_mut" = nst(has_mut & path_active),
        "PATH-_mut" = nst(has_mut & path_inactive)
    )

    names(re) = sub("PATH", pathway, names(re))
    if (length(genes) == 1)
        names(re) = sub("GENE", genes, names(re))
    else
        names(re) = sub("GENE", pathway, names(re))
    re
}

#' Wilcox test for difference in drug response between two conditions
#'
#' @param mydf  A data.frame with 
#' @param c1    The reference condition
#' @param c2    The sample condition
#' @return      A list with the reference and sample condition,
#'              p.value and fold change of median
wilcox = function(mydf, c1, c2) {
    mydf = na.omit(filter(mydf, ind %in% c(c1, c2)))
    medians = mydf %>%
        group_by(ind) %>%
        summarize(median = median(resp)) %$%
        median

    list(
        ref = c1,
        sample = c2,
        p.value = wilcox.test(resp ~ ind, data=mydf)$p.value,
        median_folds = round(10^(abs(medians[2] - medians[1])))
    )
}

#' Create a data.frame of drug responses
#'
#' @param pathway   The pathway scores
#' @param mutation  Character vector of mutated genes
#' @param drug      Drug name for which to get response IC50
#' @return          A data.frame with COSMIC IDs, subset, drug resp, etc.
create_df = function(pathway, mutation, drug) {
    strat = stratify(pathway, mutation) %>%
        stack() %>%
        transmute(cosmic=values, subset=ind)

    df$assemble(
        tissue = tissues,
        mut = apply(mut[, mutation, drop=FALSE], 1, any),
        score = scores[,pathway],
        resp = Ys[,drug]
    ) %>%
        add_rownames("cosmic") %>%
        inner_join(strat, by="cosmic") %>%
        na.omit()
}

#' Plot of the linear fit
linear_fit = function(mydf) {
    mydf %>%
        ggplot(aes(x = score, y = resp, color = tissue)) +
        geom_point(na.rm = TRUE) +
        geom_point(shape = 1, colour = 'black', na.rm = TRUE) +
        stat_smooth(method="lm", se=FALSE) +
        xlab("Pathway score") +
        ylab("Drug response") + tt
}

#' Returns a summary stat df for different contrasts
#'
#' @param start    A stratification returned by `stratify()`
#' @param pathway  A character referencing the pathway
contrast_stats = function(strat, drug, pathway, gene=pathway) {
    mydf = stack(strat)
    mydf$tissue = tissues[mydf$values]
    mydf$score = scores[,pathway][mydf$values]
    mydf$resp = Ys[,drug][mydf$values]

    as.data.frame(bind_rows(
        wilcox(mydf, paste0(gene, "_wt"), paste0(gene, "_mut")),
        wilcox(mydf, paste0(pathway, "+"), paste0(pathway, "-")),
        wilcox(mydf, paste0(pathway, "+_wt"), paste0(pathway, "-_wt")),
        wilcox(mydf, paste0(pathway, "+_mut"), paste0(pathway, "-_mut"))
    ))
}
