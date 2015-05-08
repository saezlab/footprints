YAML = $(shell find new_index -mindepth 2 -name "*.yaml")
#PATHWAYS = $(shell find new_index -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)
ACCESSION = $(basename $(notdir $(YAML)))
NORMALIZED = $(ACCESSION:%=normalized/%.RData)
ZSCORES = $(YAML:new_index/%.yaml=zscores/%.RData)

.SECONDARY:
.SECONDEXPANSION:

zscores: $(ZSCORES) zscores.RData

normalize: $(NORMALIZED)

clean:
	rm -f $(ZSCORES) $(NORMALIZED) $(NORMALIZED:%.RData=%.log)

print-%:
	@echo $* = $($*)

#zscores.RData: $(ZSCORES)
zscores.RData: $(shell find zscores -type f)
	Rscript zscores_assemble.r $@ $^

zscores/%.RData: new_index/%.yaml normalized/$$(notdir %).RData zscores_compute.r
	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@

normalized/%.RData: normalize_data.r
	@rm -f $@ # required to replace annex links
	@rm -f $(@:%.RData=%.log) # required to replace annex links
	bsub -K -M 20480 -R "rusage[mem=20480]" -R "select[gpfs]" -oo $(@:%.RData=%.log) \
		Rscript $< $(@:normalized/%.RData=%) $@
