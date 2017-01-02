MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=pathways_mapped/speed_%)

include ../../Makefile.inc
$(call ext_dep,../../model,$(notdir $(MODELS)))

.PHONY: scores

all: $(SPEED_SCORES)

# use each model to create scores object
pathways_mapped/speed_%.RData: speed.r ../../model/model_%.RData
	@mkdir -p $(dir $@)
	Rscript $< $(word 2,$^) $@
