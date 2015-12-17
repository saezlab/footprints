SCRIPTS = snp_all.r snp_drivers.r cna_gistic.r
SCORES = $(wildcard ../../scores/tcga/*.RData)

include ../../Makefile.inc
$(call ext_dep,../../scores/tcga,speed_matrix.RData)

all: $(SCRIPTS:%.r=%.pdf) $(SCORES:../../scores/tcga/%.RData=plots/%.pdf)

%.pdf: %.r ../../scores/tcga/speed_matrix.RData
	Rscript $^ $@

plots/%.pdf: pancan_heatmap.r ../../scores/tcga/%.RData
	Rscript $^ $@
