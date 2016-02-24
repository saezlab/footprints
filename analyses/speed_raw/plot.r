
pdf(OUFILE, paper="a4r")

ggplot(index, aes(x=x, y=y, color=pathway)) +
    geom_point()

dev.off()
