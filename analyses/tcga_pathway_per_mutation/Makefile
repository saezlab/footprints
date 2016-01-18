SCORES = $(wildcard ../../scores/tcga/*.RData)

include ../../Makefile.inc
$(call ext_dep,../../scores/tcga,speed_matrix.RData)

all: $(SCORES:../../scores/tcga/%.RData=heatmap/%.pdf) \
     $(SCORES:../../scores/tcga/%.RData=snp_all/%.pdf) \
     $(SCORES:../../scores/tcga/%.RData=snp_drivers/%.pdf) \
     $(SCORES:../../scores/tcga/%.RData=cna_gistic/%.pdf)

%.pdf: %.r ../../scores/tcga/speed_matrix.RData
	Rscript $^ $@

heatmap/%.pdf: pancan_heatmap.r ../../scores/tcga/%.RData
	Rscript $^ $@

snp_all/%.pdf: snp_all.r ../../scores/tcga/%.RData
	Rscript $^ $@

snp_drivers/%.pdf: snp_drivers.r ../../scores/tcga/%.RData
	Rscript $^ $@

cna_gistic/%.pdf: cna_gistic.r ../../scores/tcga/%.RData
	Rscript $^ $@

%.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
