include ../../Makefile.inc

MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
#SETS = $(wildcard ../../util/genesets/*.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)
#SET_SCORES = $(SETS:../../util/genesets/%=gsea_%)

.PHONY: prereqs scores

scores: prereqs $(SPEED_SCORES) $(SET_SCORES)

prereqs:
	make -C ../../data scores.RData
	make -C ../../model $(notdir $(MODELS))

print-%: ; @echo $* = $($*)

# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	$(call bsub_log,1024) Rscript $^ $@
