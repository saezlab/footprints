SCRIPTS = $(wildcard *.r)

all: $(SCRIPTS:%.r=%.pdf)

%.pdf: %.r ../../scores/tcga/speed_matrix.RData
	Rscript $^ $@
