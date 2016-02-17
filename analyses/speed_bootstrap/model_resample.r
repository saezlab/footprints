library(dplyr)
b = import('base')
io = import('io')
model = import('../../model/model_matrix')

ZDATA = commandArgs(TRUE)[1] %or% '../../data/zscores.RData'
OUTFILE = commandArgs(TRUE)[2] %or% 'resample1.RData'

# start from the calculated zscores.RData in ../../data
zdata = io$load(ZDATA)

# resample the index and z objects
index = zdata$index %>%
    group_by(pathway) %>%
    sample_frac(replace=TRUE) %>%
    ungroup()

zscores = zdata$zscores[,index$id]
index$id = 1:nrow(index)
colnames(zscores) = 1:ncol(zscores)

# calculate the model the same way we do with speed_matrix
model = model$zscore2model(zdata=list(index=index, zscores=zscores))

# save model resulting from resample in model/ dir
save(model, file=OUTFILE)
