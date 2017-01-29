MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
SETS = $(wildcard ../../util/genesets/mapped/*.RData)
SPEED_SCORES = $(MODEL_SCRIPTS:../../model/model_%.r=speed_%.RData)
SET_SCORES = $(SETS:../../util/genesets/mapped/%=gsea_%)
GSVA_SCORES = $(SETS:../../util/genesets/mapped/%=gsva_%)

include ../../Makefile.inc
$(call ext_dep,../../data,expr.RData)
$(call ext_dep,../../model,$(notdir $(MODEL_SCRIPTS)))
$(call ext_dep,../../util/genesets/mapped,$(notdir $(SETS)))


all: mapped

mapped: $(SPEED_SCORES) $(SET_SCORES) $(GSVA_SCORES) pathifier.RData spia.RData


# use each model to create scores object
speed_%.RData: speed.r ../../model/model_%.r ../../data/expr.RData ../../data/zscores.RData
	$(call bsub_log,40960) Rscript $^ $@

gsea_%.RData: set_enrichment.r ../../util/genesets/mapped/%.RData ../../data/expr.RData
	$(call bsub_log,4096) Rscript $^ $@

gsva_%.RData: gsva.r ../../util/genesets/mapped/%.RData ../../data/expr.RData
	$(call bsub_log,4096) Rscript $^ $@

pathifier.RData: pathifier.r ../../util/genesets/mapped/reactome.RData  ../../data/expr.RData
	$(call bsub_log,4096) Rscript $^ $@

spia.RData: spia.r ../../data/expr.RData
	$(call bsub_log,4096) Rscript $^ $@
