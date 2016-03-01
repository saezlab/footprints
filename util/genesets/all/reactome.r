b = import('base')
ar = import('array')
er = import('../enrichr')

INFILE = commandArgs(TRUE)[1] %or% "../ReactomePathways.gmt"
OUTFILE = commandArgs(TRUE)[2] %or% "reactome.RData"

pathways = er$parse_gmt(INFILE)
 
save(pathways, file=OUTFILE)
