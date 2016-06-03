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

#' Volcano plots and table of top associations
#'
#' @param fid  File ID
#' @param n    How many associations to list
volcano = function(fid, n=15) {
    fp = io$file_path('assocs_mapped', fid, ext=".RData")
    assocs = vp$load_fun(fp)$assocs.pan %>%
        arrange(p.value) %>%
        select(-std.error)
    print(vp$plot_pancan(assocs, text.size=2))
    print(head(as.data.frame(assocs), n), digits=3)
}

#' Volcano plots and table of top associations
#'
#' @param fid  File ID
#' @param n    How many associations to list
volcano_tissue = function(fid, n=15) {
    fp = io$file_path('assocs_mapped', fid, ext=".RData")
    assocs = vp$load_fun(fp)$assocs.tissue %>%
        filter(Ysub == "clinical") %>%
        mutate(Ys = ifelse(nchar(Yf) <= 13, Yf, paste0(substr(Yf, 1, 12), ">"))) %>%
        arrange(p.value) %>%
        select(-std.error, -Ysub, -Yf) %>%
        select(Ys, everything())
    print(vp$plot_pancan(assocs, p=0.1, text.size=2, base.size=5))
    print(head(as.data.frame(assocs), n), digits=3)
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
    mydf = na.omit(filter(mydf, subset %in% c(c1, c2)))
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
contrast_stats = function(strat, drug, pathway, gene=pathway) {
    mydf = stack(strat) %>%
        transmute(cosmic=values, subset=ind)
    mydf$tissue = tissues[mydf$cosmic]
    mydf$score = scores[,pathway][mydf$cosmic]
    mydf$resp = Ys[,drug][mydf$cosmic]

    as.data.frame(bind_rows(
        wilcox(mydf, paste0(gene, "_wt"), paste0(gene, "_mut")),
        wilcox(mydf, paste0(pathway, "+"), paste0(pathway, "-")),
        wilcox(mydf, paste0(pathway, "+_wt"), paste0(pathway, "-_wt")),
        wilcox(mydf, paste0(pathway, "+_mut"), paste0(pathway, "-_mut"))
    ))
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

#' Sub-stratifies the mut/wt sub-population using the (higher) pathway score
#'
#' @param path   Character string for pathway
#' @param mut    Character vector of genes whose mutations count in
#' @param drug   Character string of the drug
#' @param strat  Which sub-population should be stratified - "mut" (default) or "wt"
#' @return       A ggplot2 object with violin/boxplots
cmp_mut_path = function(path, mut, drug, strat="mut") {
    if (strat == "mut")
        levels = c("GENE_all", "GENE_wt", "GENE_mut", "PATH+_mut", "PATH-_mut", "PATH+", "PATH0", "PATH-")
    else if (strat == "wt")
        levels = c("GENE_all", "GENE_mut", "GENE_wt", "PATH+_wt", "PATH-_wt", "PATH+", "PATH0", "PATH-")
    else
        stop("need 'mut' or 'wt' for stratification")

    mydf = util$create_df(path, mut, drug) %>%
        mutate(subset = factor(subset, levels=levels)) %>%
        na.omit() # remove those subsets were we don't assign a level

    # rename GENE and PATH to what we actually look at
    levels(mydf$subset) = sub("PATH", path, levels(mydf$subset))
    levels(mydf$subset) = sub("GENE", ifelse(length(mut)==1, mut, path), levels(mydf$subset))

    ggplot(mydf, aes(x=subset, y=resp)) +
        geom_violin() +
        geom_boxplot(width=.5, outlier.shape=NA) +
        theme_bw() +
        xlab("") +
        ylab("Drug response [log uM]") #+
#        theme(axis.text.x = element_text(angle=45, hjust=1)) # makes them unevenly high
}
