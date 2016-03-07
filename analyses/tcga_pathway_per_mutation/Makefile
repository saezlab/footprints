SCORES = $(shell find ../../scores/tcga -name "*.RData")

HEATMAPS = $(SCORES:../../scores/tcga/pathways_%.RData=heatmap/%.pdf)
DRIVERS = $(SCORES:../../scores/tcga/pathways_%.RData=drivers/%.pdf)
SNP_ALL = $(SCORES:../../scores/tcga/pathways_%.RData=snp_all/%.pdf)
CNA = $(SCORES:../../scores/tcga/pathways_%.RData=cna/%.pdf)

include ../../Makefile.inc
#$(call ext_dep,../../scores/tcga,speed_matrix.RData)

all: $(HEATMAPS) $(DRIVERS) $(SNP_ALL) $(CNA)

%.pdf: %.r ../../scores/tcga/speed_matrix.RData
	Rscript $^ $@

heatmap/%.pdf: pancan_heatmap.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@

snp_all/%.pdf: snp_all.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@

drivers/%.pdf: snp_drivers.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@

cna/%.pdf: cna_gistic.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@

#%.pdf: %.Rmd
#	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
