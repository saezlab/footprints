library(dplyr)
library(grid)
library(gridExtra)
b = import('base')
io = import('io')
ar = import('array')
plt = import('plot')

venn = function() {
    # want: a list for GO/Reactome/SPEED, pathway, gene
    go = io$load(module_file('../../util/genesets/mapped/go.RData'))
    reactome = io$load(module_file('../../util/genesets/mapped/reactome.RData'))
    speed = io$load(module_file('../../model/model_matrix.RData'))$assocs %>%
        group_by(pathway) %>%
        top_n(100, -p.value) %>%
        ungroup() %>%
        unstack(gene ~ pathway) %>%
        ar$split(along=2, drop=TRUE)

    # do: Venn diagram for the overlap for all pathways
    pathways = b$list$transpose(list('Gene Ontology'=go, Reactome=reactome,
                                     'Perturbation-response'=speed))
    pathways = pathways[gtools::mixedsort(names(pathways))]

    pdf("/dev/null")

    #par(mfrow=c(3,4))
    plots = lapply(seq_along(pathways), function(i)
        grid.arrange(gTree(children = plt$vennDiagramFromList(
                            pathways[[i]], categories=c("","",""))),
                     top = textGrob(names(pathways)[i],
                     gp = gpar(fontsize=30,font=8))))

    dev.off()

    tg = grobTree(textGrob("Gene Ontology", y=0.65, vjust=0,
                           gp = gpar(fontsize=25, face=3, col="blue")),
                  textGrob("Reactome", y=0.5, vjust=0,
                           gp = gpar(fontsize=25, face=2, col="red")),
                  textGrob("Perturbation-response", y=0.35, vjust=0,
                           gp = gpar(fontsize=25, face=3, col="green")),
                  cl="titlegrob")

    plots = c(plots, list(tg))
    do.call(grid.arrange, c(plots, list(ncol=4)))
}

if (is.null(module_name())) {
    pdf("geneset_overlap_venn.pdf", width=19, height=15)
    on.exit(dev.off)

    venn()
}
