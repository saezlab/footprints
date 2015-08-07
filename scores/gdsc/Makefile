MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SETS = $(wildcard ../../util/genesets/*.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)
SET_SCORES = $(SETS:../../util/genesets/%=%)

.PHONY: prereqs scores

scores: prereqs $(SPEED_SCORES) $(SET_SCORES)

prereqs:
	make -C ../../model $(notdir $(MODELS))
	make -C ../../util/genesets $(notdir $(SETS))

print-%: ; @echo $* = $($*)

# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	Rscript $< $(word 2,$^) $@

%.RData: set_enrichment.r ../../util/genesets/%.RData ../../data/zscores.RData
	bsub -K -M 1024 -R "rusage[mem=1024]" -R "select[panfs_nobackup_research]" -oo $(@:%.RData=%.log) \
		Rscript $< $(word 2,$^) $@
