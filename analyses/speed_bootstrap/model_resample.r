b = import('base')
io = import('io')
hpc = import('hpc')

calc_resample = function(zdata, seed) {
    library(dplyr)
    set.seed(1059371 + seed)

    # resample the index and z objects
    index = zdata$index %>%
        group_by(pathway) %>%
        sample_frac(replace=TRUE) %>%
        ungroup()

    zscores = zdata$zscores[,index$id]
    index$id = 1:nrow(index)
    colnames(zscores) = 1:ncol(zscores)

    # calculate the model the same way we do with speed_matrix
    model = import('../../model/model_matrix')
    model$zscore2model(zdata=list(index=index, zscores=zscores))
}

ZDATA = commandArgs(TRUE)[1] %or% '../../data/zscores.RData'
OUTFILE = commandArgs(TRUE)[2] %or% 'model_resample.RData'
RESAMPLE = 1000

# start from the calculated zscores.RData in ../../data
zdata = io$load(ZDATA)
result = hpc$Q(calc_resample, seed=1:RESAMPLE, const=list(zdata=zdata),
               n_jobs=RESAMPLE, memory=2048)
#result = calc_resample(zdata, 0)

# save model resulting from resample in model/ dir
save(result, file=OUTFILE)
