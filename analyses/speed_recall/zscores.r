library(ggplot2)
io = import('io')
ar = import('array')

pdf("zscores.pdf", paper="a4r", width=26, height=20)
on.exit(dev.off)

# unfiltered z-scores
z = io$load('../../data/zscores.RData')
meanz = ar$map(abs(z$scores), along=1, function(x) {
    x = x[!is.na(x)]
    sum(x) / length(x)
})

# draw a histogram for all z-scores
hist(meanz, 100)

# mclust plot
library(mclust)
mcplot = Mclust(meanz)
plot(mcplot, "BIC")
plot(mcplot, "classification")
plot(mcplot, "uncertainty")
plot(mcplot, "density")

# boxplot of z-scores for each pathway
z$index$meanz = meanz
ggplot(z$index, aes(x=pathway, y=meanz)) +
    geom_boxplot() +
    geom_point() +
    xlab("pathway") +
    ylab("sum |z|") +
    ggtitle("Expression impact for different pathways")
