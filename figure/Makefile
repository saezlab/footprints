FIGURES = 1_overview.svg 2_recall.svg 3_response_expression.svg 4_drug_response.svg 5_patient_survival.svg
PDF = $(FIGURES:%.svg=%.pdf)

#.INTERMEDIATE: $(PDF)

Figures.pdf: $(PDF)
	pdftk $(PDF) cat output $@

%.pdf: %.svg
	inkscape $^ --export-pdf=$@
