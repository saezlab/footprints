library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
plt = import('plot')

# load data
scores = io$load('../../scores/gdsc/speed_matrix.RData')
Ys = gdsc$drug_response('IC50s') # or AUC
tissues = gdsc$tissues(minN=5)
ar$intersect(scores, tissues, Ys, along=1)

do_plot = function(tissue, drug, path) {
    df = data.frame(
        drug = Ys[tissues==tissue, drug],
        path = scores[tissues==tissue, path]
    )
    print(summary(lm(drug ~ path, data=df)))
    ggplot(df, aes(x=path, y=drug)) +
        geom_smooth(size=5, method=stats::lm, se=F, na.rm=T, alpha=0.1) +
        geom_point(pch=21, size=5, colour="black", alpha=1, na.rm=T) +
        theme_bw() +
        xlab(paste(path, "pathway score")) +
        ylab(paste("log IC50", drug))
}

p1 = do_plot("BRCA", "Trametinib", "MAPK")
p2 = do_plot("SCLC", "PD-0325901", "TNFa") + xlim(-1.7,0.6)
p3 = do_plot("SCLC", "FK866", "TNFa")

pcol = arrangeGrob(p1, p2, p3, ncol=1, nrow=3)

pdf("fits.pdf", width=5, height=15)
plot(pcol)
dev.off()
