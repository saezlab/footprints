library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
icgc = import('data/icgc')

OUTFILE = commandArgs(TRUE)[1] %or% "corrected_expr.h5"

# load TCGA, GDSC data for given tissue
tissues = gdsc$tissues(minN=10)
gdsc_expr = gdsc$basal_expression()

avail = icgc$available(clinical=TRUE, rna_seq=TRUE, map_to="icgc_specimen_id")
clinical = icgc$clinical(avail) %>%
    filter(grepl("Primary", specimen_type)) %>%
    select(icgc_specimen_id, project_code, tissue) %>%
    mutate(tissue = ifelse(tissue %in% c("COAD","READ"), "COREAD", tissue)) %>%
    filter(tissue %in% tissues)
tcga_expr = icgc$rna_seq(clinical$icgc_specimen_id, voom=TRUE)
tcga_expr = tcga_expr[,clinical$icgc_specimen_id] #TODO: duplicates, ordering @icgc module

tissues = tissues[tissues %in% clinical$tissue &
                  names(tissues) %in% colnames(gdsc_expr)]
gdsc_expr = gdsc_expr[,names(tissues)] #TODO: tissue filtering on expr @module

# perform batch correction
genes = intersect(rownames(tcga_expr), rownames(gdsc_expr))
expr = na.omit(cbind(tcga_expr[genes,], gdsc_expr[genes,]))
batches = c(clinical$project_code, rep("cell_line", ncol(gdsc_expr)))
covar = c(clinical$tissue, tissues)
#corrected = st$batch$combat(expr, batches, covar) #FIXME:
mat = model.matrix(~as.factor(covar))
corrected = sva::ComBat(expr, batch=batches, mod=mat, par.prior=T)

# save object with matrix for TCGA, GDSC
io$h5save(list(expr=corrected, tissue=covar), file=OUTFILE)
