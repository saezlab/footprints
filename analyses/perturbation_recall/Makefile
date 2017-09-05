SPEED = $(wildcard ../../scores/speed/*.RData)
LINCS = $(wildcard ../../scores/lincs/*.RData)
HEATMAPS = $(SPEED:../../scores/speed/%.RData=heatmap_%.pdf)
METHODS = $(shell shyaml get-values analysis_set < ../../config/methods.yaml)
METHODS2 = $(filter-out paradigm,$(METHODS)) #TODO
SCORES = $(METHODS2:%=../../scores/speed/%.RData)

plot: $(HEATMAPS) roc_speed.pdf roc_lincs.pdf speed_downsample_plot.pdf

heatmap_%.pdf: ../../scores/speed/%.RData heatmap.r
	Rscript heatmap.r $< $@

roc_speed.pdf: roc_speed.r roc_util.r sigs_gsva.RData sigs_zscores.RData $(SCORES) ../../scores/speed/speed_matrix.RData
	Rscript $< $@

roc_lincs.pdf: roc_lincs.r $(LINCS) ../../scores/lincs/speed_matrix.RData
	Rscript $< $@

sigs_gsva.RData: sigs_gsva.r ../../data/expr.RData ../../scores/speed/gsva.r
	Rscript $< $@

sigs_zscores.RData: sigs_zscores.r ../../data/expr.RData ../../data/zscores.RData ../../scores/speed/speed.r
	Rscript $< $@

speed_downsample_auc.RData: speed_downsample_auc.r ../../model/model_matrix.r ../../data/expr.RData ../../data/zscores.RData
	Rscript $^ $@

speed_downsample_plot.pdf: speed_downsample_plot.r speed_downsample_auc.RData
	Rscript $^ $@
