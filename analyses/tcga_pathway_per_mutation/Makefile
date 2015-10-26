all: cna.pdf snp.pdf snp_nodriverfilter.pdf

%.pdf: %.r ../../scores/tcga/speed_matrix.RData
	Rscript $^ $@
