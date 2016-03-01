b = import('base')
ar = import('array')
er = import('../enrichr')

INFILE = commandArgs(TRUE)[1] %or% "../Enrichr/src/main/resources/BioCarta.gmt"
OUTFILE = commandArgs(TRUE)[2] %or% "biocarta.RData"

pathways = er$parse_gmt(INFILE)
 
save(pathways, file=OUTFILE)
