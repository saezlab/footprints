include ../Makefile.inc
$(call ext_dep,../data,scores.RData)
	
MODEL_SCRIPTS = $(wildcard model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)

.PHONY: models

models: $(MODELS)

# use each model script to create a model file
model_%.RData: model_%.r ../data/scores.RData
	Rscript $< $(word 3,$^) $@
