MODEL_SCRIPTS = $(wildcard model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SCORES = $(MODELS:model_%=scores_%)
PLOTS = $(SCORES:scores_%.RData=%.pdf)

.PHONY: models scores plots

plots: $(PLOTS)

scores: $(SCORES)

models: $(MODELS)

# use each model script to create a model file
model_%.RData: model_%.r
	Rscript $< $@

# use each model to create scores object
scores_%.RData: model2scores.r model_%.RData
	Rscript $< $(word 2,$^) $@

# use each scores object to create plots
%.pdf: scores_%.RData
	Rscript $< $@
