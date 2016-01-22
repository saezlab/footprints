b = import('base')
io = import('io')
ar = import('array')
tcga = import('data/tcga')

TISSUE = commandArgs(TRUE)[1] %or% "COAD"
OUTFILE = commandArgs(TRUE)[2] %or% "paradigm_COAD_expr.txt"

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

io$write_table(expr, file=OUTFILE, sep="\t")
