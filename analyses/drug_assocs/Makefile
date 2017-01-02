#FIXME: this only works if scores target was run before
GDSC_SCORES = $(shell find ../../scores/gdsc -name "*.RData")
GDSC_ASSOCS = $(GDSC_SCORES:../../scores/gdsc/pathways_%.RData=assocs_%.RData)

CCLE_SCORES = $(shell find ../../scores/ccle -name "*.RData")
CCLE_ASSOCS = $(CCLE_SCORES:../../scores/ccle/pathways_%.RData=assocs_%_ccle.RData)

ASSOCS = $(GDSC_ASSOCS) $(CCLE_ASSOCS)
PLOTS = $(ASSOCS:assocs_%.RData=plots_%.pdf)

include ../../Makefile.inc
$(call ext_dep,../../scores/gdsc,$(notdir $(GDSC_SCORES)))
$(call ext_dep,../../scores/ccle,$(notdir $(CCLE_SCORES)))


all: mapped $(PLOTS)

mapped: $(filter plots_mapped/%,$(PLOTS)) plots_all/mutation.pdf drug_ranges.pdf


# use each model script to create a model file
$(PLOTS): plots_%.pdf: plot.r assocs_%.RData
	Rscript $^ $@

# only difference here is memory
$(filter assocs_mapped/%,$(GDSC_ASSOCS)): assocs_%.RData: gdsc_assocs.r ../../scores/gdsc/pathways_%.RData
	$(call bsub_log,4096) Rscript $^ $@

$(filter assocs_all/%,$(GDSC_ASSOCS)): assocs_%.RData: gdsc_assocs.r ../../scores/gdsc/pathways_%.RData
	$(call bsub_log,110000) Rscript $^ $@


$(filter assocs_mapped/%,$(CCLE_ASSOCS)): assocs_%_ccle.RData: ccle_assocs.r ../../scores/ccle/pathways_%.RData
	$(call bsub_log,4096) Rscript $^ $@

drug_ranges.pdf: drug_ranges.r
	Rscript $^ $@
