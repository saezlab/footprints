all: assocs_naive.pdf assocs_cov=cells.pdf assocs_4-8h.pdf

.SECONDARY:

models_%.RData: models_%.r
	bsub -K -M 10240 -oo models_$*.log -R "rusage[mem=10240]" -R "select[gpfs]" Rscript $< $@

assocs_%.pdf: models_%.RData assocs.r
	bsub -K -M 10240 -oo assocs_$*.log -R "rusage[mem=10240]" -R "select[gpfs]" Rscript assocs.r $< $@
