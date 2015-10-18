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

speed = io$load("speed_matrix.RData")

pdf(OUTFILE, paper="a4r", width=26, height=20)
print(assocfile2plot("speed_matrix.RData") + ggtitle("response genes"))
#print(assocfile2plot("gsea_reactome.RData") + ggtitle("reactome expr"))

#TODO: write this: take named vectors and sort in df
fits = df$assemble(
    FH535 = speed$Yf[,"FH535"],
    VEGF = speed$scores[,"VEGF"],
    muts = gdsc$mutated_genes(intogen=TRUE, tissue="SKCM") %>%
        ar$map(along=2, function(x) paste(names(x)[x], collapse=",")),
    FLT1 = gdsc$basal_expression()["FLT1",] # VEGFR1
)
fits$FH535 = as.numeric(fits$FH535)
fits$VEGF = as.numeric(fits$VEGF)
fits$FLT1 = as.numeric(fits$FLT1)

ggplot(fits, aes(x=VEGF, y=FH535, label=muts, fill=FLT1)) +
    geom_smooth(method=stats::lm, se=FALSE, na.rm=TRUE, alpha=0.1) +
    geom_point(pch=21, size=3, colour="black") #+
# there is not much info in mutations, FLT1 expression
#    geom_text(aes(size=0.2))
# colour by vegf expression?
dev.off()
