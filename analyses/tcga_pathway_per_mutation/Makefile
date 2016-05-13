TCGA_SCORES = $(shell find ../../scores/tcga -name "*.RData")

DRIVER_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_driver_%.RData)
CNA_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_cna_%.RData)

DRIVER_PLOTS = $(DRIVER_ASSOCS:assocs_%.RData=plot_%.pdf)
CNA_PLOTS = $(CNA_ASSOCS:assocs_%.RData=plot_%.pdf)

#DRIVER_PLOTS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=plot_driver_%.pdf)
#SNP_PLOTS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=plot_snp_%.pdf)
#CNA_PLOTS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=plot_cna_%.pdf)

HEATMAPS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=heatmap/%.pdf)

include ../../Makefile.inc
#$(call ext_dep,../../scores/tcga,speed_matrix.RData)

all: plots

# use this to no compute the _all associations, takes a bit long (and not main paper)
FILTER_OUT = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))

assocs: $(call FILTER_OUT,_all,$(ASSOCS))

$(DRIVER_ASSOCS): assocs_driver_%.RData: assocs_drivers.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	$(call bsub_log,4096) Rscript $^ $@

$(CNA_ASSOCS): assocs_cna_%.RData: assocs_cna_gistic.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	$(call bsub_log,4096) Rscript $^ $@


plots: $(call FILTER_OUT,_all,$(DRIVER_PLOTS) $(SNP_PLOTS) $(CNA_PLOTS))

$(DRIVER_PLOTS) $(CNA_PLOTS): plot_%.pdf: plot.r assocs_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@


#all: $(HEATMAPS) $(DRIVERS) $(SNP_ALL) $(CNA)
#
#%.pdf: %.r ../../scores/tcga/speed_matrix.RData
#	Rscript $^ $@
#
#heatmap/%.pdf: pancan_heatmap.r ../../scores/tcga/pathways_%.RData
#	@mkdir -p $(dir $@)
#	Rscript $^ $@
#
#snp_all/%.pdf: snp_all.r ../../scores/tcga/pathways_%.RData
#	@mkdir -p $(dir $@)
#	Rscript $^ $@
#
#drivers/%.pdf: snp_drivers.r ../../scores/tcga/pathways_%.RData
#	@mkdir -p $(dir $@)
#	Rscript $^ $@
#
#cna/%.pdf: cna_gistic.r ../../scores/tcga/pathways_%.RData
#	@mkdir -p $(dir $@)
#	Rscript $^ $@
#
#%.pdf: %.Rmd
#	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
