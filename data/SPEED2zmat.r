#!/usr/bin/env Rscript
library(modules)
gn = import('general')
bj = import('hpc/BatchJobsWrapper')

expr2zmat = function(efile, verbose=F) {
    require(modules)
    gn = import('general')
    expr = gn$loadContents(efile)
    index = gn$read.table("index.txt")
    index = index[index$GSM %in% colnames(expr) &
                  index$GPL %in% c('GPL96','GPL570','GPL571','GPL6244'),]

    subs = unique(index$SUBSET)
    scores = matrix(NA, nrow=nrow(expr), ncol=length(subs), 
                    dimnames=list(rownames(expr), subs))
    discard = c()

    for (s in subs) {
        if (verbose)
            print(s)

        negGSM = index$GSM[index$SUBSET==s & index$EFFECT=="control"]
        posGSM = index$GSM[index$SUBSET==s & index$EFFECT=="activating"]

        if (length(negGSM)>=1 && length(posGSM)>=1) {
            meanNeg = rowMeans(as.matrix(expr[,negGSM]))
            meanPos = rowMeans(as.matrix(expr[,posGSM]))
            if (length(negGSM) == 1)
                stdev = apply(expr[,c(negGSM,posGSM)], 1, sd)
            else
                stdev = apply(expr[,negGSM], 1, sd)
            
            zval = (meanPos-meanNeg) / predict(loess(stdev~meanNeg), meanPos)
#           zval = (meanPos-meanNeg) / stdev # what does loess do?
            print(head(zval))
            scores[names(zval),s] = zval
        } else
            discard = c(discard, s)
    }

    scores[is.na(scores)] = 0 # for 0-5 genes the model fails, for whatever reason
    scores[,!colnames(scores) %in% discard] # this is a lot worse for svd
}

if (is.null(module_name())) {
    # load array data
    #scores = gn$data('Raw-SPEED/zval')$scores[,subs]

    # read matrix files
    efiles = list.files(pattern="SPEED2mat_[a-zA-Z0-9_]+.RData")
    names(efiles) = gn$grepo("f?rma_[a-zA-Z0-9]+", efiles)

    #if (opt$svd) { # TODO: this needs to go into SPEED2mat
    #    for (g in unique(index$GSE)) {
    #        gsms = index$GSM[index$GSE == g]
    #        expr[,gsms] = svd(expr[,gsms])$u
    #    }
    #}

    # compute z vectors and save
    zscores = bj$Q(expr2zmat, efile=efiles, get=T)
    save(zscores, file="SPEED2zmats.RData")
}

