all: tissue_cna.pdf tissue_snp.pdf pancan_cna.pdf pancan_snp.pdf

%.pdf: %.r ../../scores/tcga/speed_matrix.RData
	Rscript $^ $@
