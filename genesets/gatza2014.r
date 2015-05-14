b = import('base')
ar = import('array')

df = read.csv('ng.3073-S2.csv', row.names=NULL, check.names=F)
lists = list(gatza=lapply(ar$split(as.matrix(df), 2), b$omit$empty))

save(lists, file = commandArgs(TRUE)[1] %or% 'gatza.RData')
