include ../../Makefile.inc

$(call ext_dep,../../scores/tcga,speed_linear.RData)
$(call ext_dep,../../scores/merge,speed_linear.RData)

all: tissue_act.pdf tact_merge.pdf

tissue_act.pdf: tissue_act.r ../../scores/tcga/speed_linear.RData
	Rscript $^ $@

tact_merge.pdf: tact_merge.r ../../scores/merge/speed_linear.RData
	Rscript $^ $@
