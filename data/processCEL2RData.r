#!/usr/bin/env R
library(modules)
gn = import('general')
ma = import('microarray')

DIR = "." #commandArgs(TRUE)[1]
NORM = commandArgs(TRUE)[2]
OUTFILE = commandArgs(TRUE)[3]
INFILES = commandArgs(TRUE)[c(-1,-2,-3)]

expr = ma$CELsToExpression(directory=DIR, method=NORM, files=INFILES)
save(expr, file=OUTFILE)

