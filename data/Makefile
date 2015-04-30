YAML = $(shell find new_index -mindepth 2 -name "*.yaml")
ACCESSION = $(basename $(notdir $(YAML)))
NORMALIZED = $(ACCESSION:%=normalized/%.RData)

normalize: $(NORMALIZED)

print-%:
	@echo $* = $($*)

normalized/%.RData: normalize_data.r
	bsub -K -M 20480 -R "rusage[mem=20480]" -R "select[gpfs]" -oo $(@:%.RData=%.log) \
		Rscript $< $(@:normalized/%.RData=%) $@
