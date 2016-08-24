library(xlsx)
library(dplyr)
b = import('base')
io = import('io')

model = io$load('../model/model_matrix.RData')$model
model = model[rowSums(model != 0) > 0,]

write.xlsx(model, file="PRG_model.xlsx", append=FALSE)
