MODELS = $(wildcard ../../scores/merge/*.RData)
PLOTS = $(MODELS:../../scores/merge/%.RData=%.pdf)

all: $(PLOTS)

%.pdf: ../../scores/merge/%.RData
	Rscript drug_assocs.r $< $@
