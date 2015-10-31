include ../../Makefile.inc

all: variant_impact.pdf

variant_impact.pdf: variant_impact.r gene_variants.RData
	Rscript $^ $@

gene_variants.RData: gene_variants.r ../tcga_pathway_per_mutation/mutations_annotated_pathwayactivities_v3_mikeformat.txt
	Rscript $^ $@
