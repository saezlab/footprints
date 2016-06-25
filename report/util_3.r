library(magrittr)
library(dplyr)
b = import('base')
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

#' Loads drug associations with mapped pathways given method id
#'
#' @param fname  Method ID
#' @return       Associations data.frame
load_fun = function(fname) {
    io$file_path('../analyses/drug_assocs/assocs_mapped', fname, ext=".RData") %>%
        io$load() %$%
        assocs.pan %>%
        filter(adj.p < 0.1) %>%
        arrange(adj.p) %>%
        mutate(num = 1:nrow(.),
               method = fname)
}

#' Returns names of a logical vector where the elements are TRUE
nst = function(x) names(x[x])

#' Generates a stratification structure using pathway and mutations
#'
#' For this, we always take the top and bottom quartiles of the pathway
#' score, either among all cell lines or the mutated (/wt) subset
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
    path_mut_null = has_mut & path > quantile(path[has_mut])[2] & path < quantile(path[has_mut])[2]
    path_wt_null = !has_mut & path > quantile(path[!has_mut])[2] & path < quantile(path[!has_mut])[2]

    re = list(
        # stratification by mutation
        "GENE_all" = nst(!is.na(path)),
        "GENE_wt" = nst(!has_mut),
        "GENE_mut" = nst(has_mut),

        # stratification by SPEED score
        "PATH+" = nst(path_active),
        "PATH-" = nst(path_inactive),
        "PATH0" = nst(path_null),

        # stratification by SPEED score, wild-type subset
        "PATH+_wt" = nst(!has_mut & path > quantile(path[!has_mut])[4]),
        "PATH-_wt" = nst(!has_mut & path < quantile(path[!has_mut])[2]),
        "PATH0_wt" = nst(path_wt_null),

        # stratification by SPEED score, mutated subset
        "PATH+_mut" = nst(has_mut & path > quantile(path[has_mut])[4]),
        "PATH-_mut" = nst(has_mut & path < quantile(path[has_mut])[2]),
        "PATH0_mut" = nst(path_mut_null)
    )
}

#' Wilcox test for difference in drug response between two conditions
#'
#' @param mydf  A data.frame with 
#' @param c1    The reference condition
#' @param c2    The sample condition
#' @return      A list with the reference and sample condition,
#'              p.value and fold change of median
wilcox = function(mydf, c1, c2) {
    mydf = filter(mydf, subset %in% c(c1, c2)) %>%
        select(cosmic, resp, subset) %>%
        na.omit()

    stopifnot(sum(duplicated(mydf$cosmic)) == 0)

    medians = mydf %>%
        group_by(subset) %>%
        summarize(median = median(resp)) %$%
        median

    list(
        ref = c1,
        sample = c2,
        p.value = wilcox.test(resp ~ subset, data=mydf)$p.value,
        median_folds = round(10^(abs(medians[2] - medians[1]))),
        n = nrow(mydf)
    )
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
contrast_stats = function(mydf) {
    gene = as.character(na.omit(unique(b$grep("^([a-zA-Z0-9]+)_wt$", mydf$subset))))
    pathway = as.character(na.omit(unique(b$grep("^([a-zA-Z0-9]+)\\+", mydf$subset))))
    as.data.frame(bind_rows(
        wilcox(mydf, paste0(gene, "_wt"), paste0(gene, "_mut")),
        wilcox(mydf, paste0(pathway, "+"), paste0(pathway, "-")),
        wilcox(mydf, paste0(pathway, "+_wt"), paste0(pathway, "-_wt")),
        wilcox(mydf, paste0(pathway, "+_mut"), paste0(pathway, "-_mut"))
    ))
}

#' Sub-stratifies the mut/wt sub-population using the (higher) pathway score
#'
#' @param path      Character string for pathway
#' @param mutation  Character vector of genes whose mutations count in
#' @param drug      Character string of the drug
#' @param strat     Which sub-pop should be stratified - "mut" (default) or "wt"
#' @return          A ggplot2 object with violin/boxplots
cmp_mut_path = function(path, mutation, drug, strat="mut") {
    if (strat == "mut")
        levels = c("GENE_all", "GENE_wt", "GENE_mut", "PATH+_mut", "PATH-_mut", "PATH+", "PATH0", "PATH-")
    else if (strat == "wt")
        levels = c("GENE_all", "GENE_mut", "GENE_wt", "PATH+_wt", "PATH-_wt", "PATH+", "PATH0", "PATH-")
    else
        stop("need 'mut' or 'wt' for stratification")

    strat = stratify(path, mutation) %>%
        stack() %>%
        transmute(cosmic=values, subset=ind)

    mydf = df$assemble(
        tissue = tissues,
        mut = apply(mut[, mutation, drop=FALSE], 1, any),
        score = scores[,path],
        resp = Ys[,drug]
    ) %>%
        add_rownames("cosmic") %>%
        inner_join(strat, by="cosmic")

    # rename GENE and PATH to what we actually look at
    mydf$label = mydf$subset
    mydf$subset = as.character(mydf$subset)
    mydf$subset = sub("PATH", path, mydf$subset)
    mydf$subset = sub("GENE", ifelse(length(mutation)==1, mutation, path), mydf$subset)
    mydf$label = factor(mydf$label, levels=levels)
    levels(mydf$label) = sub("PATH", path, levels(mydf$label))
    levels(mydf$label) = sub("GENE", ifelse(length(mutation)==1, mutation, path), levels(mydf$label))
    mydf
}

plot_mut_path = function(mydf) {
    ggplot(na.omit(mydf), aes(x=label, y=resp)) +
        geom_violin() +
        geom_boxplot(width=.5, outlier.shape=NA) +
        theme_bw() +
        xlab("") +
        ylab("Drug response [log uM]") #+
#        theme(axis.text.x = element_text(angle=45, hjust=1)) # makes them unevenly high
}
