TCGA_SCORES = $(shell find ../../scores/tcga -name "*.RData")

DRIVER_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_driver_%.RData)
SNP_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_snp_%.RData)
CNA_ASSOCS = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=assocs_cna_%.RData)

DRIVER_PLOTS = $(SCORES:../../scores/tcga/pathways_%.RData=drivers/%.pdf)
SNP_PLOTS = $(SCORES:../../scores/tcga/pathways_%.RData=snp_all/%.pdf)
CNA_PLOTS = $(SCORES:../../scores/tcga/pathways_%.RData=cna/%.pdf)

HEATMAPS = $(SCORES:../../scores/tcga/pathways_%.RData=heatmap/%.pdf)

include ../../Makefile.inc
#$(call ext_dep,../../scores/tcga,speed_matrix.RData)

all: assocs


assocs: $(DRIVER_ASSOCS) $(SNP_ASSOCS) $(CNA_ASSOCS)

$(DRIVER_ASSOCS): assocs_driver_%.RData: assocs_drivers.r ../../scores/tcga/pathways_%.RData
	Rscript $^ $@

$(SNP_ASSOCS): assocs_snp_%.RData: assocs_snp.r ../../scores/tcga/pathways_%.RData
	Rscript $^ $@

$(CNA_ASSOCS): assocs_cna_%.RData: assocs_cna_gistic.r ../../scores/tcga/pathways_%.RData
	Rscript $^ $@


plots: $(DRIVER_PLOTS) $(SNP_PLOTS) $(CNA_PLOTS)


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
