MODEL_SCRIPTS = $(wildcard ../../model/model_*.RData)
MODELS = $(MODEL_SCRIPTS:%.RData=%.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=pathways_mapped/speed_%)

SETS_MAPPED = $(wildcard ../../util/genesets/mapped/*.RData)
SETS_ALL = $(wildcard ../../util/genesets/all/*.RData)
SET_SCORES_MAPPED = $(SETS_MAPPED:../../util/genesets/mapped/%=pathways_mapped/gsea_%)
SET_SCORES_ALL = $(SETS_ALL:../../util/genesets/all/%=pathways_all/gsea_%)
GSVA_SCORES_MAPPED = $(SETS_MAPPED:../../util/genesets/mapped/%=pathways_mapped/gsva_%)
GSVA_SCORES_ALL = $(SETS_ALL:../../util/genesets/all/%=pathways_all/gsva_%)

include ../../Makefile.inc
$(call ext_dep,../../model,$(notdir $(MODELS)))
$(call ext_dep,../../util/genesets/mapped,$(notdir $(SETS_MAPPED)))
$(call ext_dep,../../util/genesets/all,$(notdir $(SETS_ALL)))


all: mapped $(SPEED_SCORES) $(SET_SCORES_ALL) $(GSVA_SCORES_ALL) pathways_all/pathifier.RData pathways_all/spia.RData pathways_mapped/paradigm.RData

mapped: $(SET_SCORES_MAPPED) $(GSVA_SCORES_MAPPED) pathways_mapped/pathifier.RData pathways_mapped/spia.RData pathways_mapped/paradigm.RData


# mapped pathways
pathways_mapped/speed_%.RData: speed.r ../../model/model_%.RData ../../config/tcga.yaml
	Rscript $< $(word 2,$^) $@

pathways_mapped/gsea_%.RData: set_enrichment.r ../../util/genesets/mapped/%.RData ../../config/tcga.yaml
	$(call bsub_log,30720) Rscript $< $(word 2,$^) $@

pathways_mapped/gsva_%.RData: gsva.r ../../util/genesets/mapped/%.RData ../../config/tcga.yaml
	$(call bsub_log,30720) Rscript $< $(word 2,$^) $@

pathways_mapped/spia.RData: spia.r ../../config/tcga.yaml
	$(call bsub_log,30720) Rscript $< $@ TRUE

pathways_mapped/pathifier.RData: pathifier.r ../../util/genesets/mapped/reactome.RData ../../config/tcga.yaml
	$(call bsub_log,30720) Rscript $< $(word 2,$^) $@

pathways_mapped/paradigm.RData: paradigm.r ../../config/tcga.yaml
	Rscript $< $@

# all pathways
pathways_all/gsea_%.RData: set_enrichment.r ../../util/genesets/all/%.RData ../../config/tcga.yaml
	$(call bsub_log,40960) Rscript $< $(word 2,$^) $@

pathways_all/gsva_%.RData: gsva.r ../../util/genesets/all/%.RData ../../config/tcga.yaml
	$(call bsub_log,40960) Rscript $< $(word 2,$^) $@

pathways_all/spia.RData: spia.r ../../config/tcga.yaml
	$(call bsub_log,61140) Rscript $< $@ FALSE

pathways_all/pathifier.RData: pathifier.r ../../util/genesets/all/reactome.RData ../../config/tcga.yaml
	$(call bsub_log,61140) Rscript $< $(word 2,$^) $@

#paradigm.RData: paradigm.r
#	Rscript $^ $@
