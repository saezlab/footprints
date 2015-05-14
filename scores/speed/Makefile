MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SCORES = $(MODELS:../../model/model_%=%)

.PHONY: models scores

scores: $(SCORES)

print-%: ; @echo $* = $($*)

../../model/model_%.RData:
	make -C $(basename $@) $@

# use each model to create scores object
%.RData: speed.r ../../model/model_%.RData
	Rscript $< $(word 2,$^) $@
