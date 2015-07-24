#GENESETS = $(wildcard ../../genesets/*.RData)
#SETSCORES = $(GENESETS:../../genesets/%=%)
MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)

.PHONY: prereqs

all: prereqs $(SPEED_SCORES) pathifier.RData spia.RData #rppa.RData $(SETSCORES) $(SPEED_SCORES)

prereqs:
	make -C ../../model $(notdir $(wildcard ../../model/*.RData))
	make -C ../../expr_cluster corrected_expr.h5

rppa.RData: rppa.r
	Rscript $<

pathifier.RData: pathifier.r ../../genesets/reactome.RData ../../expr_cluster/corrected_expr.h5
	Rscript $^ $@

spia.RData: spia.r ../../genesets/reactome.RData ../../expr_cluster/corrected_expr.h5
	Rscript $^ $@

speed_%.RData: speed.r ../../model/model_%.RData ../../expr_cluster/corrected_expr.h5
	bsub -K -R "select[panfs_nobackup_research]" -oo $(@:%.RData=%.log) \
		Rscript $< $(word 2,$^) $(word 3,$^) $@

#%.RData: set_enrichment.r ../../genesets/%.RData ../../data/zscores.RData
#	@rm -f $@ # annex link overwrite
#	bsub -K -M 10240 -R "rusage[mem=10240]" -R "select[panfs_nobackup_research]" -oo $(@:%.RData=%.log) \
#		Rscript $< $(word 2,$^) $@
