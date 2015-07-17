library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
tcga = import('data/tcga')

OUTFILE = commandArgs(TRUE)[1] %or% "corrected_expr.h5"

# load GDSC data
gdsc_expr = t(gdsc$basal_expression())
gdsc_tissues = gdsc$tissues(minN=10)
ar$intersect(gdsc_expr, gdsc_tissues, along=1)
gdsc_expr = t(gdsc_expr)

# load TCGA data
tcga_expr = ar$stack(lapply(tcga$tissues(), tcga$rna_seq), along=2)
tcga_index = tcga$barcode2index(colnames(tcga_expr)) %>%
    transmute(id = Bio.ID,
              study = Study.Abbreviation,
              type = Sample.Definition) %>%
    mutate(tissue = ifelse(study %in% c("COAD","READ"), "COREAD", study),
           covar = paste0(tissue, ifelse(type=="Solid Tissue Normal", "_N", ""))) %>%
    filter(tissue %in% unique(gdsc_tissues))

# intersect tissues and subset matrices/index
tcga_expr = tcga_expr[,tcga_index$id]

# define covariates and batches
batches = c(tcga_index$study, rep("cell_line", ncol(gdsc_expr)))
covar = unname(c(tcga_index$covar, gdsc_tissues))

# perform batch correction
ar$intersect(tcga_expr, gdsc_expr, along=1)
expr = cbind(tcga_expr, gdsc_expr)
mat = model.matrix(~as.factor(covar))
set.seed(1562681)
corrected = sva::ComBat(expr, batch=batches, mod=mat, par.prior=TRUE)

# save object with matrix for TCGA, GDSC
io$h5save(list(expr=t(corrected), tissue=covar), file=OUTFILE)
