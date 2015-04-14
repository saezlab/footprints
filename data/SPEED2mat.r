#!/usr/bin/env Rscript
library(optparse)
library(stringr)
library(modules)
gn = import('general')
ar = import('array')
batch = import('batch')

optionList = list(
    make_option("--norm", help="rma, frma or soft", default="rma"),
    make_option("--method", help="combat, dwd, or none", default="combat"),
    make_option("--batch", help="gpl or gse", default="gpl")
)
opt = parse_args(OptionParser(option_list=optionList))

# load SPEED index file
index = gn$read.table("index.txt")

if (opt$norm == "soft") {
    expr = gn$data("SPEED-Data/SPEED2_soft")

    valid = intersect(index$GSM[index$EFFECT %in% c('control','activating') &
                                index$GPL %in% c('GPL570','GPL571','GPL96','GPL6244')], colnames(expr))

    expr = expr[,valid]
    idsub = unique(index[c('GPL','GSE','GSM')])
    if (opt$batch == "gpl") {
        batches = unlist(sapply(valid, function(x) idsub[idsub$GSM==x,'GPL']))
    } else if (opt$batch == "gse") {
        batches = unlist(sapply(valid, function(x) idsub[idsub$GSM==x,'GSE']))
    }
} else {
    # load array data
    datadir = "../DATA/RAW/Transcriptomic/SPEED/CELS"
    regex = paste0("GSE[0-9]+_GPL(570|571|96|6244)_", opt$norm, "\\.RData")
    arrays = gn$loadFilesByRegex(regex, datadir)

    # process matrices
    arrays = gn$delete.NULLs(lapply(arrays, function(x) {
        x = x[,colnames(x) %in% index$GSM]
        effectsInGSM = index$EFFECT[index$GSM %in% colnames(x)]
        if (all(c("control","activating") %in% effectsInGSM))
            x
        else
            NULL
    }))
    expr = ar$stack(arrays, along=2)
}

expr = na.omit(expr) #TODO: option for imputation

if (opt$norm != "soft") { #FIXME: this is horrible, find common solution for all
#TODO: get this from index, not arrays structure
    GSMboth = intersect(index$GSM, colnames(expr))
    if (opt$batch == "gpl") {
        fact = gn$grepo("GPL[0-9]+", names(arrays))
        batches = unlist(lapply(seq_along(arrays), function(i) rep(fact[i], ncol(arrays[[i]]))))
    } else if (opt$batch == "gse")
        batches = unlist(lapply(seq_along(arrays), function(i) rep(i, ncol(arrays[[i]]))))
}

expr = batch$merge(opt$method, expr, batch=batches, covariate=NULL)

# save resulting object
save(expr, file=paste0("SPEED2mat_", opt$norm, "_", paste0(opt$method,opt$batch), ".RData"))

