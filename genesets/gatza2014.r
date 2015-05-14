b = import('base')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% 'ng.3073-S2.csv'

df = read.csv('ng.3073-S2.csv', row.names=NULL, check.names=FALSE)
lists = lapply(ar$split(as.matrix(df), 2, drop=TRUE), b$omit$empty)

save(lists, file = commandArgs(TRUE)[2] %or% 'gatza.RData')
