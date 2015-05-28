YAML = $(shell find new_index -mindepth 2 -name "*.yaml")
ACCESSION = $(basename $(notdir $(YAML)))
NORMALIZED = $(ACCESSION:%=normalized/%.RData)
SCORES = $(YAML:new_index/%.yaml=scores/%.RData)

.SECONDARY:
.SECONDEXPANSION:
.PHONY: scores normalize clean

scores: $(SCORES) zscores.RData dscores.RData

normalize: $(NORMALIZED)

clean:
	rm -f $(ZSCORES) $(NORMALIZED) $(NORMALIZED:%.RData=%.log)

print-%:
	@echo $* = $($*)

%scores.RData: scores_assemble.r $(SCORES)
	@rm -f $@ # required to replace annex links
	Rscript $< $* $@ $(filter-out $<,$^)

scores/%.RData: new_index/%.yaml normalized/$$(notdir %).RData scores_compute.r
	@mkdir -p $(shell dirname $@)
	Rscript scores_compute.r $(word 1,$^) $(word 2,$^) $@

normalized/%.RData: normalize_data.r
	@mkdir -p $(shell dirname $@)
	bsub -K -M 20480 -R "rusage[mem=20480]" -oo $(@:%.RData=%.log) \
		Rscript $< $(@:normalized/%.RData=%) $@
