TXT = $(wildcard *.txt)

all: index.RData expr.RData

index.RData: $(TXT)
	Rscript create_index.r

expr.RData: index.RData
	Rscript create_expr_matrix.r
