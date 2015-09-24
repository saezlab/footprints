MODELS = $(wildcard ../../scores/merge/*.RData)
PLOTS = $(MODELS:../../scores/merge/%.RData=%.pdf)

all: $(PLOTS)

%.pdf: ../../scores/merge/%.RData drug_assocs.r
	Rscript drug_assocs.r $< $@
