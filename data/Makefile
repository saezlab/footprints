YAML = $(shell find new_index -mindepth 2 -name "*.yaml")
#PATHWAYS = $(shell find new_index -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)
ACCESSION = $(basename $(notdir $(YAML)))
NORMALIZED = $(ACCESSION:%=normalized/%.RData)
ZSCORES = $(YAML:new_index/%.yaml=zscores/%.RData)

.SECONDARY:

zscores: zscores.RData

normalize: $(NORMALIZED)

clean:
	rm -f $(ZSCORES) $(NORMALIZED) $(NORMALIZED:%.RData=%.log)

print-%:
	@echo $* = $($*)

#zscores.RData: $(ZSCORES)
zscores.RData: $(shell find zscores -type f)
	Rscript zscores_assemble.r $@ $^

#zscores/$(1)/%.RData: new_index/$(1)/%.yaml normalized/%.RData zscores_compute.r
#	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@

zscores/EGFR/%.RData: new_index/EGFR/%.yaml normalized/%.RData zscores_compute.r
	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@

zscores/MAPK/%.RData: new_index/MAPK/%.yaml normalized/%.RData zscores_compute.r
	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@

zscores/PI3K/%.RData: new_index/PI3K/%.yaml normalized/%.RData zscores_compute.r
	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@

zscores/p53/%.RData: new_index/p53/%.yaml normalized/%.RData zscores_compute.r
	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@

zscores/Trail/%.RData: new_index/Trail/%.yaml normalized/%.RData zscores_compute.r
	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@

zscores/H2O2/%.RData: new_index/H2O2/%.yaml normalized/%.RData zscores_compute.r
	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@

zscores/Hypoxia/%.RData: new_index/Hypoxia/%.yaml normalized/%.RData zscores_compute.r
	Rscript zscores_compute.r $(word 1,$^) $(word 2,$^) $@


normalized/%.RData: normalize_data.r
	@rm -f $@ # required to replace annex links
	@rm -f $(@:%.RData=%.log) # required to replace annex links
	bsub -K -M 20480 -R "rusage[mem=20480]" -R "select[gpfs]" -oo $(@:%.RData=%.log) \
		Rscript $< $(@:normalized/%.RData=%) $@
