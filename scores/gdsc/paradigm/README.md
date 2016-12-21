PARADIGM scores for GDSC data
=============================

Scripts in this directory process individual GDSC samples using the PARADIGM
tool.

The workflow is defined by make, so just process by typing:

```bash
make -j500
```

This will create the directories `expr` and `path` (explained below). File
system usage should be moderate, because preprocessing only parallelises
tissues and the file sizes for the samples are small. However, do adjust this
if the file system can not handle it (like `-j30`). Pathway scores will only be
computed once all expression data has been processed.

To remove all generated files type:

```bash
make clean
```

### Splitting GDSC gene expression in separate text files

First, the scripts `preprocess_expr.r` takes each tissue gene expression,
discards genes that are not in the PARADIGM network (`proteins.txt`,
derived from `SuperPathway.txt`), and
writes all remaining expression values in the directory
`expr/<tissue>/<barcode>.txt`

Expression values are converted into ranks as suggested by the [PARADIGM
documentation](https://sysbio.soe.ucsc.edu/paradigm/ParadigmDocumentation.pdf).

This way, we can process each sample separately using PARADIGM because
otherwise the tools runs for too long.

### Calculating pathway scores

For each of the sample files in the `expr/` directory, we use PARADIGM to
process all samples to pathway scores in the `path/` directory. This has the
same structure `<tissue>/<barcode>.txt`.

PARADIGM settings used are saved in the file `config.txt`.
