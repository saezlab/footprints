SCRIPTS = $(wildcard *.r)

include ../../Makefile.inc
$(call ext_dep,../../scores/tcga,speed_matrix.RData)

all: $(SCRIPTS:%.r=%.pdf)

%.pdf: %.r ../../scores/tcga/speed_matrix.RData
	Rscript $^ $@
