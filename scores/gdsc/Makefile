MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SETS = $(wildcard ../../genesets/*.RData)
SCORES = $(MODELS:../../model/model_%=speed_%) $(SETS:../../genesets/%=%)

.PHONY: models scores

scores: $(SCORES)

print-%: ; @echo $* = $($*)

# call other makefiles for dependencies
../../model/model_%.RData:
	make -C $(shell dirname $@) $@

../../genesets/%.RData:
	make -C $(shell dirname $@) $@

# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	Rscript $< $(word 2,$^) $@

%.RData: set_enrichment.r ../../genesets/%.RData
	Rscript $< $(word 2,$^) $@
