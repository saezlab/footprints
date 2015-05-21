GENESETS = $(wildcard ../../genesets/*.RData)
SETSCORES = $(GENESETS:../../genesets/%=%)

.PHONY: $(GENESETS)

all: speed_linear.RData speed_matrix.RData

../../genesets/%.RData:
	make -C $(shell dirname $@) $@

prereqs:
	make -C ../../model $(notdir $(wildcard ../../model/*.RData))
	make -C ../../genesets $(notdir $(GENESETS))

speed_%.RData: speed.r prereqs ../../model/model_%.RData
	@rm -f $@ # annex link overwrite
	bsub -K -M 10240 -R "rusage[mem=10240]" -R "select[gpfs]" -oo $(@:%.RData=%.log) \
		Rscript $< $(word 3,$^) $@

%.RData: set_enrichment.r prereqs ../../genesets/%.RData
	@rm -f $@ # annex link overwrite
	bsub -K -M 10240 -R "rusage[mem=10240]" -R "select[gpfs]" -oo $(@:%.RData=%.log) \
		Rscript $< $(word 3,$^) $@
