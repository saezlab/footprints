include ../../Makefile.inc

MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SETS = $(wildcard ../../util/genesets/*.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)
SET_SCORES = $(SETS:../../util/genesets/%=gsea_%)

.PHONY: prereqs scores

scores: prereqs $(SPEED_SCORES) $(SET_SCORES) pathifier.RData spia.RData

prereqs:
	make -C ../../data scores.RData
	make -C ../../model $(notdir $(MODELS))
	make -C ../../util/genesets $(notdir $(SETS))

print-%: ; @echo $* = $($*)

# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	Rscript $^ $@

gsea_%.RData: set_enrichment.r ../../util/genesets/%.RData ../../data/scores.RData
	Rscript $^ $@

pathifier.RData: pathifier.r ../../util/genesets/reactome.RData  ../../data/expr.RData
	Rscript $^ $@

spia.RData: spia.r ../../data/expr.RData
	Rscript $^ $@
