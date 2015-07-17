`corrected_expr.h5`

file where overlapping gene expression data between the GDSC and the TCGA is batch-merged using ComBat

hdf5 file has 2 objects:

* `expr` - expression matrix, samples x genes
* `tissues` - which tissues each sample is (`_N` suffix for normals)
