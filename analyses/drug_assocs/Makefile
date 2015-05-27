# this only works if scores target was run before
SCORES = $(wildcard ../../scores/gdsc/*.RData)
PLOTS = $(SCORES:../../scores/gdsc/%.RData=%.pdf)

.PHONY: prereqs

all: prereqs $(PLOTS)

prereqs:
	make -C ../../scores/gdsc $(notdir $(SCORES))

print-%: ; @echo $* = $($*)

# use each model script to create a model file
%.pdf: drug_assocs.r ../../scores/gdsc/%.RData
	Rscript $< $(word 2,$^) $@
