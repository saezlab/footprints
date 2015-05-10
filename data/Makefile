YAML = $(shell find new_index -mindepth 2 -name "*.yaml")
ACCESSION = $(basename $(notdir $(YAML)))
NORMALIZED = $(ACCESSION:%=normalized/%.RData)
ZSCORES = $(YAML:new_index/%.yaml=zscores/%.RData)

.SECONDARY:
.SECONDEXPANSION:
.PHONY: zscores normalize clean

zscores: $(ZSCORES) zscores.RData

normalize: $(NORMALIZED)

dirs:
	@mkdir -p normalized $(shell dirname $(ZSCORES))

clean:
	rm -f $(ZSCORES) $(NORMALIZED) $(NORMALIZED:%.RData=%.log)

print-%:
	@echo $* = $($*)

zscores.RData: $(ZSCORES)
	@rm -f $@ # required to replace annex links
	Rscript zscores_assemble.r $@ $^

zscores/%.RData: new_index/%.yaml normalized/$$(notdir %).RData zscores_compute.r
	@mkdir -p $(shell dirname $@)
	@rm -f $@ # required to replace annex links
	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@

normalized/%.RData: normalize_data.r
	@mkdir -p $(shell dirname $@)
	@rm -f $@ # required to replace annex links
	@rm -f $(@:%.RData=%.log) # required to replace annex links
	bsub -K -M 20480 -R "rusage[mem=20480]" -R "select[gpfs]" -oo $(@:%.RData=%.log) \
		Rscript $< $(@:normalized/%.RData=%) $@
