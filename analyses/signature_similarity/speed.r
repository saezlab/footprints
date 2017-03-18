library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% '../../scores/speed/speed_matrix.RData'
OUTFILE = commandArgs(TRUE)[2] %or% "/dev/null"

data = io$load(INFILE)
scores = t(data$scores * data$index$sign)

save(scores, file=OUTFILE)
