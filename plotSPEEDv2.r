#!/usr/bin/env Rscript
library(modules)
gn = import('general')

load('SPEEDv2.ro')

index = gn$read.table("index.txt")

fcIdx = gn$reIndexDf(index, 'GSM', colnames(fc))
plot(tsne(t(na.omit(fc))),
     col=as.factor(fcIdx$PATHWAY),
     pch=as.integer(as.factor(fcIdx$GSE)))

deltaIdx = gn$reIndexDf(index, 'GSM', colnames(delta))
plot(tsne(t(na.omit(delta))),
     col=as.factor(deltaIdx$PATHWAY),
     pch=as.integer(as.factor(deltaIdx$GSE)))

#TODO: distribute weights for gse-tissue-pathway combinations (but: how?)
#TODO: compare the quality of models we can build with delta/fc raw/weighted svm/psvm
#  (basically, just run CVd mlr

library(mlr)

