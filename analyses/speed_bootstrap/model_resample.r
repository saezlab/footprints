b = import('base')
io = import('io')
hpc = import('hpc')

calc_resample = function(zdata, zfun, seed_offset) {
    library(dplyr)
    set.seed(1059371 + seed_offset)

    # resample the index and z objects
    index = zdata$index %>%
        group_by(pathway) %>%
        sample_frac(replace=TRUE) %>%
        ungroup()

    zscores = zdata$zscores[,index$id]
    index$id = 1:nrow(index)
    colnames(zscores) = 1:ncol(zscores)

    # calculate the model the same way we do with speed_matrix
    re = zfun(zdata=list(index=index, zscores=zscores))
    re$model[rowSums(model != 0) != 0,]
}

ZDATA = commandArgs(TRUE)[1] %or% '../../data/zscores.RData'
MODULE = commandArgs(TRUE)[2] %or% '../../model/model_matrix.r'
OUTFILE = commandArgs(TRUE)[3] %or% 'model_resample.RData'

# start from the calculated zscores.RData in ../../data
zdata = io$load(ZDATA)
zfun = import_(sub("\\.r$", "", MODULE))$zscore2model

result = hpc$Q(calc_resample, seed_offset = 1:1000,
               const = list(zdata=zdata, zfun = zfun),
               n_jobs = 1000, memory = 4096)

# save model resulting from resample in model/ dir
save(result, file=OUTFILE)
