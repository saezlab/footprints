all: expr.RData plot.pdf

expr.RData: expr.r
	Rscript expr.r

plot.pdf: do.r expr.RData
	Rscript do.r
