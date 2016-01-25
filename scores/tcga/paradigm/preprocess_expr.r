b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')

TISSUE = commandArgs(TRUE)[1] %or% "COAD"

ranks = function(x) {
    re = sort(x)
    re = setNames(1:length(re) / length(re), names(re))
    re[names(x)]
}

proteins = io$read_table(module_file("proteins.txt"))$V1

expr = tcga$rna_seq(TISSUE) %>%
    ar$map(along=2, ranks) %>%
    t()

expr = expr[,intersect(colnames(expr), proteins)]

expr_df = cbind(id = rownames(expr),
                as.data.frame(expr))
rownames(expr_df) = 1:nrow(expr_df)

for (sample in 1:nrow(expr_df))
    io$write_table(expr_df[sample,,drop=FALSE],
                   file = io$file_path(TISSUE, expr_df$id[sample], ext=".txt"),
                   sep = "\t")
