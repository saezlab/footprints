library(dplyr)
.ar = import('array')
.gdsc = import('data/gdsc')

.tissues = .gdsc$tissues()
.Ys = .gdsc$drug_response('IC50s')
.ar$intersect(.tissues, .Ys, along=1)

.min_conc = .gdsc$drug$conc('min', colnames(.Ys), log=TRUE)
.max_conc = .gdsc$drug$conc('max', colnames(.Ys), log=TRUE)

#' @param ...    Parameters passed to plt$volcano
drug_tissue_volcano = function() {
}

#' @param drug      Name of the drug
#' @param stratify  Either a character vector with tissues to highlight, or
#'                  A nested list with each item a (names, for the x tick) list
#'                  of COSMIC IDs to include in stratification
#' @param min_n     Minumum number of drug reponse points to include in plot;
#'                  stratify'd points will always be plotted
drug_range_box = function(drug, stratify=NULL, min_n=5) {
    mydf = data.frame(tissue=.tissues, cosmic = names(.tissues), drug=.Ys[,drug]) %>%
        na.omit() %>%
        group_by(tissue) %>%
        filter(n() >= min_n) %>%
        ungroup() %>%
        mutate(fill = "other")

    for (i in seq_along(stratify)) {
        iobj = stratify[i]
        # get name of list item for stratified highlights
        # or element for simple highlights
        iname = names(stratify)[i] %or% iobj
        if (is.list(iname))
            stop("stratify needs to be named list or character vector")

        # highlight the requested tissue
        mydf[mydf$tissue == iname,'fill'] = iname

        # this is the stratification per tissue
        if (is.list(iobj)) {
            for (j in seq_along(iobj)) {
                jobj = iobj[[j]]
                jname = names(iobj)[j]

                # this is the stratification within tissue
                for (k in seq_along(jobj)) {
                    kobj = jobj[[k]]
                    kname = names(jobj)[k]

                    # add a stratified version of the drug response to mydf
                    mydf = mydf %>%
                        filter(cosmic %in% kobj) %>%
                        mutate(tissue = kname,
                               fill = jname) %>%
                        bind_rows(mydf)
                }
            }
        }
    }

    minc = .min_conc[drug]
    maxc = .max_conc[drug]
    rect = data.frame(xmin=-Inf, xmax=Inf, ymin=minc, ymax=maxc)

    ggplot(mydf, aes(x=reorder(tissue, drug,
            FUN=function(x) median(x, na.rm=TRUE)), y=drug, fill=fill)) +
        scale_x_discrete() +
        geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                  fill="plum2", alpha=0.1, inherit.aes=FALSE) +
        geom_boxplot(na.rm=TRUE) +
        guides(fill=FALSE) +
        xlab("Cancer type") +
        ylab("IC50 [log uM]") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        geom_boxplot(na.rm=TRUE) +
        geom_abline(intercept=minc, slope=0, linetype="dotted") +
        geom_abline(intercept=maxc, slope=0, linetype="dotted") +
        ggtitle(paste("Drug response ranges per tissue for", drug))
}

drug_fit = function() {
}

if (is.null(module_name())) {
# plot some example plot
}
