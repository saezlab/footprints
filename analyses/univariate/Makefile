# this only works if scores target was run before
SCORES = $(wildcard ../../scores/speed/*.RData)
PLOTS = $(SCORES:../../scores/speed/%.RData=%.pdf)

.PHONY: models scores

plots: $(PLOTS)

print-%: ; @echo $* = $($*)

# use each model script to create a model file
%.pdf: drug_assocs.r ../../scores/speed/%.RData
	Rscript $< $(word 2,$^) $@
