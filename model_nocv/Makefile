include ../Makefile.inc
	
MODEL_SCRIPTS = $(wildcard model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)

.PHONY: models prereqs

models: prereqs $(MODELS)

prereqs:
	make -C ../data scores.RData

# use each model script to create a model file
model_%.RData: model_%.r ../data/scores.RData
	Rscript $< $(word 3,$^) $@
