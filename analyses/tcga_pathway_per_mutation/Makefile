all: tissue.pdf pancan.pdf

%.pdf: %.r ../../scores/tcga/speed_matrix.RData
	Rscript $^ $@
