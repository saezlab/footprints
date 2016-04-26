include ../Makefile.inc
$(call ext_dep,../data,zscores.RData zscores_commonNULL.RData)

MODEL_SCRIPTS = $(wildcard model_*.r)
MODELS = $(MODEL_SCRIPTS:%.r=%.RData)

.PHONY: models

all: $(MODELS) #model_cmat.RData model_cmat_robust.RData

# use each model script to create a model file
model_%.RData: model_%.r ../data/zscores.RData
	Rscript $^ $@

#model_cmat.RData: model_matrix.r ../data/zscores_commonNULL.RData
#	Rscript $^ $@
#
#model_cmat_robust.RData: model_matrix_robust.r ../data/zscores_commonNULL.RData
#	Rscript $^ $@
