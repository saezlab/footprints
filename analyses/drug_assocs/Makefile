# this only works if scores target was run before
SCORES = $(wildcard ../../scores/gdsc/*.RData)
PLOTS = $(SCORES:../../scores/gdsc/%.RData=%.pdf)

.PHONY: models prereqs

plots: $(PLOTS)

prereqs:
	make -C ../../scores/gdsc $(notdir $(SCORES))

print-%: ; @echo $* = $($*)

# use each model script to create a model file
%.pdf: drug_assocs.r prereqs ../../scores/gdsc/%.RData
	Rscript $< $(word 3,$^) $@
