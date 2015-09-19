include ../../Makefile.inc

MODEL_SCRIPTS = $(wildcard ../../model_nocv/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SPEED_SCORES = $(MODELS:../../model_nocv/model_%=speed_%)

.PHONY: prereqs scores

scores: prereqs $(SPEED_SCORES) $(SET_SCORES)

prereqs:
	make -C ../../data scores.RData
	make -C ../../model $(notdir $(MODELS))

# use each model to create scores object
speed_%.RData: speed.r ../../model_nocv/model_%.RData
	$(call bsub_log,1024) Rscript $^ $@
