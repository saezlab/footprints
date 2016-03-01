b = import('base')
ar = import('array')
er = import('../enrichr')

INFILE = commandArgs(TRUE)[1] %or% "../Enrichr/src/main/resources/ChEA.gmt"
OUTFILE = commandArgs(TRUE)[2] %or% "chea.RData"

lists = er$parse_gmt(INFILE)
 
save(lists, file=OUTFILE)
