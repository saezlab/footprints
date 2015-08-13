include ../../Makefile.inc

#GENESETS = $(wildcard ../../util/genesets/*.RData)
#SETSCORES = $(GENESETS:../../util/genesets/%=%)
MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)

.PHONY: prereqs

all: prereqs $(SPEED_SCORES) pathifier.RData spia.RData #rppa.RData $(SETSCORES) $(SPEED_SCORES)

prereqs:
	make -C ../../model $(notdir $(wildcard ../../model/*.RData))
	make -C ../../util/expr_cluster corrected_expr.h5

rppa.RData: rppa.r
	Rscript $<

pathifier.RData: pathifier.r ../../util/genesets/reactome.RData ../../util/expr_cluster/corrected_expr.h5
	Rscript $^ $@

spia.RData: spia.r ../../util/genesets/reactome.RData ../../util/expr_cluster/corrected_expr.h5
	Rscript $^ $@

speed_%.RData: speed.r ../../model/model_%.RData ../../util/expr_cluster/corrected_expr.h5
	$(call bsub_log) Rscript $< $(word 2,$^) $(word 3,$^) $@

#%.RData: set_enrichment.r ../../util/genesets/%.RData ../../data/zscores.RData
#	@rm -f $@ # annex link overwrite
#	bsub -K -M 10240 -R "rusage[mem=10240]" -R "select[panfs_nobackup_research]" -oo $(@:%.RData=%.log) \
#		Rscript $< $(word 2,$^) $@
