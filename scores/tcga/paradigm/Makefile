# TODO: pull this in over config/yaml
TISSUES = ACC BLCA BRCA CESC CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD

all: $(TISSUES:%=paradigm_%_expr.txt)

proteins.txt:
	grep ^protein < SuperPathway.txt | grep -v abstract | cut -f2 > $@

paradigm_%_expr.txt:
	Rscript preprocess_expr.r $* $@
