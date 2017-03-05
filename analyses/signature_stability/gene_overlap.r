library(dplyr)
library(ggplot2)
io = import('io')
ar = import('array')

# ~ 10 min
nums = io$load('./resample_matrix.RData') %>%
    ar$stack(along=3, fill=0) %>%
    ar$map(along=3, function(x) sum(x!=0)) %>%
    reshape2::melt() %>%
    select(gene = Var1,
           pathway = Var2,
           n = value) %>%
    filter(n != 0)

p = ggplot(nums, aes(x=n)) +
    geom_histogram(binwidth=1) +
    facet_wrap(~pathway, scales="free")

if (is.null(module_name())) {
    pdf("gene_overlap.pdf", paper="a4r")
    print(p)
    dev.off()
}
