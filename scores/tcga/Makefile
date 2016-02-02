GENESETS = $(wildcard ../../util/genesets/*.RData)
SETSCORES = $(GENESETS:../../util/genesets/%=gsea_%)
MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)

include ../../Makefile.inc
$(call ext_dep,../../model,$(notdir $(wildcard ../../model/*.RData)))
$(call ext_dep,../../util/genesets,$(notdir $(GENESETS)))


all: pathifier.RData spia.RData paradigm.RData $(SETSCORES) $(SPEED_SCORES) # rppa.RData

rppa.RData: rppa.r
	Rscript $<

spia.RData: spia.r
	$(call bsub_log,1024) Rscript $^ $@

paradigm.RData: paradigm.r
	Rscript $^ $@

pathifier.RData: pathifier.r ../../util/genesets/reactome.RData
	$(call bsub_log,1024) Rscript $^ $@

speed_%.RData: speed.r ../../model/model_%.RData
	$(call bsub_log,10240) Rscript $^ $@

gsea_%.RData: set_enrichment.r ../../util/genesets/%.RData
	$(call bsub_log,10240) Rscript $^ $@
