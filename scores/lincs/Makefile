MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SETS = $(wildcard ../../util/genesets/mapped/*.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)
GSVA_SCORES = $(SETS:../../util/genesets/mapped/%=gsva_%)

include ../../Makefile.inc
$(call ext_dep,../../util/lincs,index.RData)
$(call ext_dep,../../model,$(notdir $(MODELS)))
$(call ext_dep,../../util/genesets/mapped,$(notdir $(SETS)))


all: mapped

mapped: $(SPEED_SCORES) $(GSVA_SCORES) pathifier.RData spia.RData


# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.RData
	$(call bsub_log,4096) Rscript $^ $@

gsea_%.RData: set_enrichment.r ../../util/genesets/mapped/%.RData ../../util/lincs/index.RData
	$(call bsub_log,4096) Rscript $^ $@

gsva_%.RData: gsva.r ../../util/genesets/mapped/%.RData ../../util/lincs/index.RData
	$(call bsub_log,4096) Rscript $^ $@

pathifier.RData: pathifier.r ../../util/genesets/mapped/reactome.RData  ../../util/lincs/index.RData
	$(call bsub_log,4096) Rscript $^ $@

spia.RData: spia.r ../../util/lincs/index.RData
	$(call bsub_log,4096) Rscript $^ $@
