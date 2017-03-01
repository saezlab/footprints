all: expr.RData plot.pdf

expr.RData: expr.r validation.yaml
	Rscript expr.r

plot.pdf: do.r expr.RData
	Rscript do.r
