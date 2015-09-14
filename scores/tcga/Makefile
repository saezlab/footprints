include ../../Makefile.inc
	
GENESETS = $(wildcard ../../util/genesets/*.RData)
SETSCORES = $(GENESETS:../../util/genesets/%=gsea_%)
MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)

.PHONY: prereqs

all: prereqs pathifier.RData spia.RData $(SETSCORES) $(SPEED_SCORES) # rppa.RData

prereqs:
	make -C ../../model $(notdir $(wildcard ../../model/*.RData))
	make -C ../../util/genesets $(notdir $(GENESETS))

rppa.RData: rppa.r
	Rscript $<

spia.RData: spia.r ../../util/genesets/reactome.RData
	$(call bsub_log,1024) Rscript $^ $@

pathifier.RData: pathifier.r ../../util/genesets/reactome.RData
	$(call bsub_log,1024) Rscript $^ $@

speed_%.RData: speed.r ../../model/model_%.RData
	$(call bsub_log,10240) Rscript $^ $@

gsea_%.RData: set_enrichment.r ../../util/genesets/%.RData ../../data/zscores.RData
	$(call bsub_log,10240) Rscript $^ $@
