library(modules)
b = import('base')
hpc = import('hpc')

OUTFILE = commandArgs(TRUE)[1] %or% "models_4-8h.RData"

trainModel = function(nfactors) {
    library(peer)
    library(modules)
    library(dplyr)
    io = import('io')
    ar = import('array')

    # load speed dscores
    dscores = io$data('SPEED-Data/SPEED2dmats')$rma_none

    # set up prior
    index = io$read_table("../SPEED-Data/zval_meta_BTOmapped.txt", header=T) %>%
        filter(effect == "activating") %>%
        filter(time %in% c("4 h","5 h","6 h","5-6 h","5.66 h","7 h","8 h")) %>%
        select(id, pathway, cells) %>%
        filter(id %in% colnames(dscores))

    dscores = dscores[,index$id]

    model=PEER() # initialize peer
    PEER_setNk(model, nfactors)#ncol(dscores)) # number of hidden factors (just large enough)
    PEER_setPhenoMean(model, dscores) # phenotype = expression data
    #PEER_setCovariates(model, as.matrix(covariates)) # optional: set covariates
    PEER_setNmax_iterations(model, 10000) # optional: set number of steps; default: 1000
    PEER_update(model) # calculate the model

    # calc drug assocs w/ hidden factors
    # optional/later: GO-validate them
    w = PEER_getW(model) # weights of different factors in different experiments
    x = PEER_getX(model)
    rownames(x) = rownames(dscores)
    rownames(w) = colnames(dscores)
    list(x=x, w=w)
}

nfactors = c(5,10,15,20,25,30,40,50,60,80,100,150)
models = hpc$Q(trainModel, nfactors=setNames(nfactors, nfactors))
save(models, file=OUTFILE)
