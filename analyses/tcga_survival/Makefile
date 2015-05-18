TCGA_SCORES = $(wildcard ../../scores/tcga/*.RData)
SURV = $(TCGA_SCORES:../../scores/tcga/%.RData=%.pdf)

all: $(SURV)

../../scores/tcga/%.RData:
	make -C $(shell dirname $@) $@

%.pdf: ../../scores/tcga/%.RData survival.r
	Rscript survival.r $< $@
