GENESETS = $(wildcard ../../genesets/*.RData)
SETSCORES = $(GENESETS:../../genesets/%=%)
MODEL_SCRIPTS = $(wildcard ../../model/model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)
SPEED_SCORES = $(MODELS:../../model/model_%=speed_%)

.PHONY: prereqs

all: $(SPEED_SCORES) $(SETSCORES) pathifier.RData spia.RData #rppa.RData $(SETSCORES) $(SPEED_SCORES)

#prereqs:
#	make -C ../../model $(notdir $(wildcard ../../model/*.RData))
#	make -C ../../expr_cluster corrected_expr.h5

#rppa.RData: rppa.r
#	Rscript $<

pathifier.RData: pathifier.r ../../genesets/reactome.RData
	Rscript $^ $@

spia.RData: spia.r ../../data/lincs_perturbation_qc/index.RData
	Rscript $^ $@

speed_%.RData: speed.r ../../model/model_%.RData ../../data/lincs_perturbation_qc/index.RData
	Rscript $^ $@

%.RData: set_enrichment.r ../../genesets/%.RData ../../data/lincs_perturbation_qc/index.RData
	Rscript $^ $@
