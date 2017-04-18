library(dplyr)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% '../../scores/speed/speed_matrix.RData'
OUTFILE = commandArgs(TRUE)[2] %or% "/dev/null"

data = io$load(INFILE)

# we already scale=1 in scores
# this is to show relative pathway activations
scores = t(data$scores * data$index$sign) %>%
    ar$map(along=2, scale)

save(scores, file=OUTFILE)
