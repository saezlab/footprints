b = import('base')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% '../ng.3073-S2.csv'
OUTFILE = commandArgs(TRUE)[2] %or% 'gatza.RData'

df = read.csv(INFILE, row.names=NULL, check.names=FALSE)
pathways = lapply(ar$split(as.matrix(df), 2, drop=TRUE), b$omit$empty)

save(pathways, file=OUTFILE)
