MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SETS = $(wildcard ../../util/genesets/*.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)
SET_SCORES = $(SETS:../../util/genesets/%=gsea_%)

include ../../Makefile.inc
$(call ext_dep,../../data,expr.RData)
$(call ext_dep,../../model,$(notdir $(MODELS)))
$(call ext_dep,../../util/genesets,$(notdir $(SETS)))


all: $(SPEED_SCORES) $(SET_SCORES) pathifier.RData spia.RData

# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	$(call bsub_log,4096) Rscript $^ $@

gsea_%.RData: set_enrichment.r ../../util/genesets/%.RData ../../data/expr.RData
	$(call bsub_log,4096) Rscript $^ $@

pathifier.RData: pathifier.r ../../util/genesets/reactome.RData  ../../data/expr.RData
	$(call bsub_log,4096) Rscript $^ $@

spia.RData: spia.r ../../data/expr.RData
	$(call bsub_log,4096) Rscript $^ $@
