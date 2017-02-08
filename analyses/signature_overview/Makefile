include ../../Makefile.inc

all: geneset_overlap.pdf geneset_overlap_venn.pdf dataset_size_bars.pdf gene_significance.pdf speed_treemap.pdf speed_go_plot.pdf

%.pdf: %.r
	Rscript $< $@

dataset_size_bars.pdf: dataset_size_bars.r dataset_size.txt
	Rscript $< $@

speed_go_piano.RData: speed_go_piano.r
	$(call bsub_log,10240) Rscript $< $@

speed_go_plot.pdf: speed_go_plot.r speed_go_piano.RData
	Rscript $< $@
