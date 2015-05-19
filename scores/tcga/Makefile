GENESETS = $(wildcard ../../genesets/*.RData)
SETSCORES = $(GENESETS:../../genesets/%=%)

.PHONY: $(GENESETS)

all: speed_linear.RData speed_matrix.RData $(SETSCORES)

../../genesets/%.RData:
	make -C $(shell dirname $@) $@

speed_%.RData: speed.r ../../model/model_%.RData
	Rscript $< $(word 2,$^) $@

%.RData: set_enrichment.r ../../genesets/%.RData
	bsub -K -M 10240 -R "rusage[mem=10240]" -R "select[gpfs]" -oo $(@:%.RData=%.log) \
	Rscript $< $(word 2,$^) $@
