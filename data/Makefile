include ../Makefile.inc

YAML = $(shell find new_index -mindepth 2 -name "*.yaml")
ACCESSION = $(basename $(notdir $(YAML)))
NORMALIZED = $(ACCESSION:%=normalized/%.RData)
SCORES = $(YAML:new_index/%.yaml=scores/%.RData)

#.SECONDARY:
#.SECONDEXPANSION:
.PHONY: scores normalize clean

scores: zscores.RData

expr: expr.RData

normalize: $(NORMALIZED)

clean:
	rm -f $(SCORES) $(NORMALIZED) $(NORMALIZED:%.RData=%.log)

expr.RData: expr.r $(NORMALIZED)
	Rscript $^

zscores.RData: zscores.r expr.RData
	Rscript $^ $@

normalized/%.RData: normalize_data.r
	@mkdir -p $(shell dirname $@)
	$(call bsub_log,20480) Rscript $< $(@:normalized/%.RData=%) $@
