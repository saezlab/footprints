TCGA_SCORES = $(shell find ../../scores/tcga -name "*.RData")

CONT_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_cont_%.RData)
DISCR_ASSOCS = #$(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_discr_%.RData)

CONT_PLOTS = $(CONT_ASSOCS:assocs_%.RData=plot_%.pdf)
DISCR_PLOTS = #$(DISCR_ASSOCS:assocs_%.RData=plot_%.pdf)

FILTER_OUT = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))
ASSOCS = 
include ../../Makefile.inc
#$(call ext_dep,../../scores/tcga,$(notdir $(CONT) $(DISCR)))

all: plots

assocs: $(call FILTER_OUT,_all,$(DISCR_ASSOCS) $(CONT_ASSOCS))

$(CONT_ASSOCS): assocs_cont_%.RData: assocs_cont.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	$(call bsub_log,4096) Rscript $^ $@

$(DISCR_ASSOCS): assocs_discr_%.RData: assocs_discr.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	$(call bsub_log,4096) Rscript $^ $@


plots: $(call FILTER_OUT,_all,$(DISCR_PLOTS) $(CONT_PLOTS))

$(CONT_PLOTS) $(DISCR_PLOTS): plot_%.pdf: plot.r assocs_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@
