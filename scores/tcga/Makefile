include ../../Makefile.inc
	
GENESETS = $(wildcard ../../util/genesets/*.RData)
SETSCORES = $(GENESETS:../../util/genesets/%=%)
MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)

.PHONY: prereqs

all: prereqs rppa.RData $(SETSCORES) $(SPEED_SCORES)

prereqs:
	make -C ../../model $(notdir $(wildcard ../../model/*.RData))
	make -C ../../util/genesets $(notdir $(GENESETS))

rppa.RData: rppa.r
	Rscript $<

speed_%.RData: speed.r ../../model/model_%.RData
	$(call bsub_log,10240) Rscript $< $(word 2,$^) $@

%.RData: set_enrichment.r ../../util/genesets/%.RData ../../data/zscores.RData
	$(call bsub_log,10240) Rscript $< $(word 2,$^) $@
