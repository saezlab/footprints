MODEL_SCRIPTS = $(wildcard model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)

.PHONY: models

models: $(MODELS)

print-%: ; @echo $* = $($*)

# use each model script to create a model file
model_%.RData: model_%.r ../data/zscores.RData
	Rscript $< $(word 2,$^) $@
