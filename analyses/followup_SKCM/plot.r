library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
df = import('data_frame')
plt = import('plot')
gdsc = import('data/gdsc')

INFILE = commandArgs(TRUE)[1] %or% "../../scores/gdsc/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "SKCM_assocs.pdf"

assocfile2plot = function(fname) {
    assocs.tissue = io$load(fname)$assocs.tissue

    # volcano plot for tissue subsets
    assocs.tissue %>%
        mutate(label = paste(Yf, scores, sep=":")) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", dir=-1) %>%
        plt$volcano(p=0.2)
}

do_fit = function(tissue, drug, path, gene) {
    fits = df$assemble(
        drug = speed$Yf[,drug],
        path = speed$scores[,path],
        muts = gdsc$mutated_genes(intogen=TRUE, tissue=tissue) %>%
            ar$map(along=2, function(x) paste(names(x)[x], collapse=",")),
        gene = gdsc$basal_expression()[gene,]
    )
    fits$drug = as.numeric(fits$drug)
    fits$path = as.numeric(fits$path)
    fits$gene = as.numeric(fits$gene)

    ggplot(fits, aes(x=path, y=drug, label=muts, fill=gene)) +
        geom_smooth(method=stats::lm, se=FALSE, na.rm=TRUE, alpha=0.1) +
        geom_point(pch=21, size=5, colour="black") +
        scale_fill_gradient(low="white", high="black",
                            limits=c(min(fits$gene, na.rm=TRUE),
                                     max(fits$gene, na.rm=TRUE)),
                            guide = guide_legend(title = gene)) +
        theme_bw() +
        ggtitle("SKCM: FH535 response and VEGF activity") +
        xlab("predicted VEGF activity") +
        ylab("log IC50 FH535 [beta-catenin]")
}

speed = io$load("speed_matrix.RData")

pdf(OUTFILE, paper="a4r", width=26, height=20)
print(assocfile2plot(INFILE))
print(do_fit("SKCM", "FH535", "VEGF", "PDGFRA"))
print(do_fit("SKCM", "FH535", "VEGF", "PDGFRB"))
dev.off()
