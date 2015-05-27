GENESETS = $(wildcard ../../genesets/*.RData)
SETSCORES = $(GENESETS:../../genesets/%=%)

.PHONY: prereqs

all: prereqs speed_linear.RData speed_matrix.RData

prereqs:
	make -C ../../model $(notdir $(wildcard ../../model/*.RData))
	make -C ../../genesets $(notdir $(GENESETS))

speed_%.RData: speed.r ../../model/model_%.RData
	@rm -f $@ # annex link overwrite
	bsub -K -M 10240 -R "rusage[mem=10240]" -R "select[gpfs]" -oo $(@:%.RData=%.log) \
		Rscript $< $(word 2,$^) $@

%.RData: set_enrichment.r ../../genesets/%.RData
	@rm -f $@ # annex link overwrite
	bsub -K -M 10240 -R "rusage[mem=10240]" -R "select[gpfs]" -oo $(@:%.RData=%.log) \
		Rscript $< $(word 2,$^) $@
