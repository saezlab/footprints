b = import('base')

INFILE = commandArgs(TRUE)[1] %or% "speed_go.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "speed_go_plot.pdf"
