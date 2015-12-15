include ../../Makefile.inc

# what if import targets could match pattern instead of
# individual files that makefile creates
# --> this could save quite a bit of var def

TCGA_SCORES = $(wildcard ../../scores/tcga/*.RData) #TODO: this only updates, does not create
$(call import_target, ../../scores/tcga, $(TCGA_SCORES))
#TCGA_SCORES = $(call import_target_pattern, ../../scores/tcga, %.RData) # this would be preferred
PLOTS = $(TCGA_SCORES:../../scores/tcga/%.RData=plots/%.pdf)

all: $(PLOTS)

assocs/%.RData: assocs.r ../../scores/tcga/%.RData gene_variants.RData
	$(call bsub_log,4096) Rscript $^ $@

plots/%.pdf: plot.r assocs/%.RData
	Rscript $^ $@

gene_variants.RData: gene_variants.r ../tcga_pathway_per_mutation/mutations_annotated_pathwayactivities_v3_mikeformat.txt
	Rscript $^ $@

%.pdf: %.Rmd util.r
	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
