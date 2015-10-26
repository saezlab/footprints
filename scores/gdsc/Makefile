include ../../Makefile.inc

MODEL_SCRIPTS = $(wildcard ../../model/model_*.RData)
MODELS = $(MODEL_SCRIPTS:%.RData=%.RData)
SETS = $(wildcard ../../util/genesets/*.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)
SET_SCORES = $(SETS:../../util/genesets/%=gsea_%)

.PHONY: prereqs scores

scores: prereqs $(SPEED_SCORES) $(SET_SCORES)

prereqs:
	make -C ../../model $(notdir $(MODELS))
	make -C ../../util/genesets $(notdir $(SETS))

# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	Rscript $< $(word 2,$^) $@

gsea_%.RData: set_enrichment.r ../../util/genesets/%.RData ../../data/expr.RData
	$(call bsub_log,1024) Rscript $< $(word 2,$^) $@
