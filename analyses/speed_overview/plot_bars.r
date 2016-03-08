do_plot = function() {   
    library(cowplot)
    library(reshape2)
    io = import('io')

    nums = read.table(module_file("numbers.txt"), header=TRUE)
    nums$speed2 = NULL
    nums = melt(nums, id="num")
    nums$num = factor(nums$num, levels=c("pathways", "experiments", "contrasts", "arrays"))

    ggplot(nums, aes(x=variable, y=value, fill=variable)) +
        geom_bar(stat="identity") +
        facet_wrap(~ num, scales="free") +
        xlab("Signature Method") +
        ylab("Number") +
        theme(axis.text.x = element_text(angle=45, hjust=1),
              legend.position = "none")
}

if (is.null(module_name())) {
    pdf("bars.pdf")

    do_plot()

    dev.off()
}
