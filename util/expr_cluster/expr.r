library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
tcga = import('data/tcga')

OUTFILE = commandArgs(TRUE)[1] %or% "corrected_expr.RData"

# specify the tissues that we want
include_tissues = import('../../config')$tcga$tissues_with_normals

# load GDSC data
gdsc_expr = gdsc$basal_expression()
gdsc_index = data_frame(id=colnames(gdsc_expr)) %>%
    mutate(study = "gdsc",
           type = "Cancer Cell Line",
           tissue = gdsc$cosmic$id2tissue(id)) %>%
    mutate(covar = tissue,
           batch = "gdsc_cell_line") %>%
    filter(tissue %in% include_tissues | tissue == "COREAD") %>%
    distinct()

# load TCGA data
tcga_expr = ar$stack(lapply(include_tissues, tcga$rna_seq), along=2)
tcga_index = tcga$barcode2index(colnames(tcga_expr)) %>%
    transmute(id = Bio.ID,
              study = Study.Abbreviation,
              type = Sample.Definition) %>%
    mutate(tissue = ifelse(study %in% c("COAD","READ","COADREAD"), "COREAD", study),
           covar = paste0(tissue, ifelse(type=="Solid Tissue Normal", "_normal", "")),
           batch = study) %>%
    distinct()

# assemble data objects
tcga_expr = tcga_expr[,tcga_index$id]
gdsc_expr = gdsc_expr[,gdsc_index$id]
ar$intersect(tcga_expr, gdsc_expr, along=1)
raw_expr = cbind(tcga_expr, gdsc_expr)
index = bind_rows(tcga_index, gdsc_index)

# make sure we don't have duplicates and ids match
stopifnot(sum(duplicated(index$id)) == 0)
stopifnot(index$id == colnames(raw_expr))

# create model matrix and perform batch correction
mat = model.matrix(~as.factor(index$covar))
set.seed(1562681)
expr = sva::ComBat(raw_expr, batch=index$batch, mod=mat, par.prior=TRUE)

# save object with matrix for TCGA, GDSC
save(expr, index, file=OUTFILE)
