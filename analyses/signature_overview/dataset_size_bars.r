library(cowplot)
library(reshape2)
io = import('io')

create_df = function() {
    nums = read.table(module_file("dataset_size.txt"), header=TRUE)
    levels = nums$measure
    nums$speed2 = NULL
    nums = melt(nums, id="measure")
    nums$measure = factor(nums$measure, levels=levels)
    nums
}

do_plot = function(nums) {   
    ggplot(nums, aes(x=variable, y=value, fill=variable)) +
        geom_bar(stat="identity") +
        facet_wrap(~ measure, scales="free") +
        xlab("Signature Method") +
        ylab("Number") +
        theme(axis.text.x = element_text(angle=45, hjust=1),
              legend.position = "none")
}

if (is.null(module_name())) {
    pdf("dataset_size_bars.pdf")

    print(do_plot(create_df()))

    dev.off()
}
