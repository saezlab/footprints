MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SETS = $(wildcard ../../genesets/*.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)
SET_SCORES = $(SETS:../../genesets/%=%)

.PHONY: prereqs scores

scores: prereqs $(SPEED_SCORES) $(SET_SCORES)

prereqs:
	make -C ../../model $(notdir $(MODELS))
	make -C ../../genesets $(notdir $(SETS))

print-%: ; @echo $* = $($*)

# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	Rscript $< $(word 2,$^) $@

%.RData: set_enrichment.r ../../genesets/%.RData
	@rm -f $@ # annex link overwrite
	bsub -K -M 1024 -R "rusage[mem=1024]" -oo $(@:%.RData=%.log) \
		Rscript $< $(word 2,$^) $@
