library(dplyr)
library(bhtsneR)
b = import('base')
io = import('io')
ar = import('array')

INFILE = commandArgs(TRUE)[1] %or% '../../scores/speed/speed_fullmat.RData'
OUTFILE = commandArgs(TRUE)[2] %or% "/dev/null"

data = io$load(INFILE)
index = data$index

scores = data$scores
scores[is.na(scores)] = 0 #FIXME:

stopifnot(rownames(scores) == index$id)
dim2 = tsne(scores)

index = index %>%
    mutate(x = dim2[,1],
           y = dim2[,2])

save(index, file=OUTFILE)
