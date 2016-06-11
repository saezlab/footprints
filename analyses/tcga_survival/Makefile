TCGA_SCORES = $(shell find ../../scores/tcga -name "*.RData")

CONT_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_cont_%.RData)
DISCR_ASSOCS = #$(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_discr_%.RData) #discrete currently failing

CONT_PLOTS = $(CONT_ASSOCS:assocs_%.RData=plot_%.pdf)
DISCR_PLOTS = #$(DISCR_ASSOCS:assocs_%.RData=plot_%.pdf)

include ../../Makefile.inc
$(call ext_dep,../../scores/tcga,$(notdir $(CONT) $(DISCR)))

MAPPED = $(call contains,_mapped,$(DISCR_PLOTS) $(CONT_PLOTS))


all: $(CONT_PLOTS) $(DISCR_PLOTS)

mapped: $(call contains,_mapped,$(DISCR_PLOTS) $(CONT_PLOTS))


$(CONT_ASSOCS): assocs_cont_%.RData: assocs_cont.r util.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	$(call bsub_log,4096) Rscript $< $(word 3,$^) $@

$(DISCR_ASSOCS): assocs_discr_%.RData: assocs_discr.r util.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	$(call bsub_log,4096) Rscript $< $(word 3,$^) $@


$(CONT_PLOTS) $(DISCR_PLOTS): plot_%.pdf: plot.r assocs_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@
