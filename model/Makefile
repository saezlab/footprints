include ../Makefile.inc
	
MODEL_SCRIPTS = $(wildcard model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)

.PHONY: models prereqs

models: prereqs $(MODELS)

prereqs:
	make -C ../data expr.RData zscores.RData

# use each model script to create a model file
model_%.RData: model_%.r ../data/zscores.RData
	Rscript $^ $@
