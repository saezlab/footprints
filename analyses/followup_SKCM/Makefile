include ../../Makefile.inc
	
ASSOCS = speed_matrix gsea_reactome
ASSOC_IN = $(ASSOCS:%=../../scores/gdsc/%.RData)
ASSOC_OUT = $(ASSOCS:%=%.RData)

all: $(ASSOC_OUT)

#all: SKCM_assocs.pdf

#SCKM_assocs.pdf: plot.r $(ASSOC_OUT)
#	Rscript plot.r

speed_matrix.RData: drug_assocs.r ../../scores/gdsc/speed_matrix.RData
	Rscript $^ $@

gsea_reactome.RData: drug_assocs.r ../../scores/gdsc/gsea_reactome.RData
	Rscript $^ $@
