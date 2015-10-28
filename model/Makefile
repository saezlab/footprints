include ../Makefile.inc
$(call ext_dep,../data,zscores.RData)

MODEL_SCRIPTS = $(wildcard model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)

.PHONY: models

models: $(MODELS)

# use each model script to create a model file
model_%.RData: model_%.r ../data/zscores.RData
	Rscript $^ $@
