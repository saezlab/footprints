MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)

include ../../Makefile.inc
$(call ext_dep,../../model,$(notdir $(MODELS)))

.PHONY: scores

scores: $(SPEED_SCORES) $(SET_SCORES)

# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	Rscript $< $(word 2,$^) $@
