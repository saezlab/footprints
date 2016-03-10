TCGA_SCORES = $(shell find ../../scores/tcga -name "*.RData")

CONT_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_cont_%.RData)
DISCR_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_discr_%.RData)

CONT_PLOTS = $(CONT_ASSOCS:assocs_%=plot_%.pdf)
DISCR_PLOTS = $(DISCR_ASSOCS:assocs_%=plot_%.pdf)


include ../../Makefile.inc
#$(call ext_dep,../../scores/tcga,$(notdir $(CONT) $(DISCR)))

all: assocs


assocs: $(CONT_ASSOCS) $(DISCR_ASSOCS)

$(CONT_ASSOCS): assocs_cont_%.RData: assocs_cont.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@

$(DISCR_ASSOCS): assocs_discr_%.RData: assocs_discr.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@


plots: $(CONT_PLOTS) $(DISCR_PLOTS)

$(CONT_PLOTS): plot_cont_%.pdf: plot_cont.r assocs_cont%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@

$(DISCR_PLOTS): plot_discr_%.pdf: plot_discr.r assocs_discr%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@
