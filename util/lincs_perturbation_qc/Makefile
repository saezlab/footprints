TXT = $(wildcard *.txt)

all: index.RData expr.RData

index.RData: create_index.r $(TXT)
	Rscript create_index.r

expr.RData: create_expr_matrix.r index.RData
	Rscript create_expr_matrix.r
