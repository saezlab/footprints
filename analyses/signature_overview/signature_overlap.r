library(dplyr)
b = import('base')
io = import('io')
plt = import('plot')

overlap = function() {
    df = io$load(module_file('../../model/model_matrix.RData'))$assocs %>%
        group_by(pathway) %>%
        top_n(100, -p.value) %>%
        ungroup() %>%
        select(gene, pathway)

    pathways = gtools::mixedsort(unique(df$pathway))

    p2num = function(p1, p2)
        length(intersect(df$gene[df$pathway==p1], df$gene[df$pathway == p2]))

    mat = matrix(NA, nrow=length(pathways), ncol=length(pathways),
        dimnames=list(rev(pathways), pathways))
    for (p1 in pathways)
        for (p2 in pathways)
            mat[p1, p2] = p2num(p1, p2)

    plt$matrix(reshape2::melt(mat), value ~ Var1 + Var2, palette="Blues") +
        geom_text(aes(label=value)) +
        coord_fixed() +
        theme(legend.position = "none") +
        xlab("") +
        ylab("")
}

if (is.null(module_name())) {
    pdf("signature_overlap.pdf")
    print(overlap())
    dev.off()
}
