include ../../Makefile.inc

all: model_resample.RData

model_resample.RData: model_resample.r ../../data/zscores.RData
	Rscript $^ $@
