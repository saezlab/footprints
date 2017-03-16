library(xlsx)
library(dplyr)
b = import('base')
io = import('io')

INFILE = commandArgs(TRUE)[1] %or% "../model/model_matrix.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "PROGENy_model.xlsx"

model = io$load(INFILE)$model
model = model[rowSums(model != 0) > 0,]

write.xlsx(model, file=OUTFILE, append=FALSE)
