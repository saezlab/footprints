TCGA_SCORES = $(shell find ../../scores/tcga -name "*.RData")

DRIVER_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_driver_%.RData)
CNA_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_cna_%.RData)

DRIVER_PLOTS = $(DRIVER_ASSOCS:assocs_%.RData=plot_%.pdf)
CNA_PLOTS = $(CNA_ASSOCS:assocs_%.RData=plot_%.pdf)


include ../../Makefile.inc
$(call ext_dep,../tcga_mutation,cna_driver_matrix.RData)


all: mapped #$(DRIVER_PLOTS) $(CNA_PLOTS)

mapped: $(call contains,_mapped,$(DRIVER_PLOTS) $(SNP_PLOTS) $(CNA_PLOTS))


$(DRIVER_ASSOCS): assocs_driver_%.RData: assocs_drivers.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	$(call bsub_log,10240) Rscript $^ $@

$(CNA_ASSOCS): assocs_cna_%.RData: assocs_cna_gistic.r ../../scores/tcga/pathways_%.RData ../tcga_mutation/cna_driver_matrix.RData
	@mkdir -p $(dir $@)
	$(call bsub_log,10240) Rscript $^ $@


$(DRIVER_PLOTS) $(CNA_PLOTS): plot_%.pdf: plot.r assocs_%.RData
	@mkdir -p $(dir $@)
	$(call bsub_log,10240) Rscript $^ $@
