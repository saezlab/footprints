MODEL_SCRIPTS = $(wildcard ../../model/model_*.RData)
MODELS = $(MODEL_SCRIPTS:%.RData=%.RData)
SETS = $(wildcard ../../util/genesets/*.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)
SET_SCORES = $(SETS:../../util/genesets/%=gsea_%)

include ../../Makefile.inc
$(call ext_dep,../../model,$(notdir $(MODELS)))
$(call ext_dep,../../util/genesets,$(notdir $(SETS)))

.PHONY: scores

all: $(SPEED_SCORES) $(SET_SCORES)

# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	Rscript $< $(word 2,$^) $@

gsea_%.RData: set_enrichment.r ../../util/genesets/%.RData
	$(call bsub_log,1024) Rscript $^ $@

spia.RData: spia.r
	Rscript $^ $@

pathifier.RData: pathifier.r ../../util/genesets/reactome.RData
	Rscript $^ $@

paradigm.RData: paradigm.r
	Rscript $^ $@
