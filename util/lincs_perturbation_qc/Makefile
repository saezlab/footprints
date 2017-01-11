TXT = $(wildcard *.txt)

all: index.RData expr.RData

index.RData: index.r $(TXT)
	Rscript index.r $@ $^

expr.RData: expr.r index.RData
	Rscript expr.r
