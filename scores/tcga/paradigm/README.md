PARADIGM scores for TCGA data
=============================

Scripts in this directory process individual TCGA samples using the PARADIGM
tool.

The workflow is defined by make, so just process by typing:

```bash
make -j500 # processing a sample takes a long time, so run 500 in parallel
```

This will create the directories `expr` and `path` (explained below).

To remove all generated files type:

```bash
make clean
```

### Splitting TCGA gene expression in separate text files

First, the scripts `preprocess_expr.r` takes each tissue gene expression,
discards genes that are not in the PARADIGM network (`proteins.txt`,
derived from `SuperPathway.txt`), and
writes all remaining expression values in the directory
`expr/<tissue>/<barcode>.txt`

This way, we can process each sample separately using PARADIGM because
otherwise the tools runs for too long.

### Calculating pathway scores

For each of the sample files in the `expr/` directory, we use PARADIGM to
process all samples to pathway scores in the `path/` directory. This has the
same structure `<tissue>/<barcode>.txt`.

PARADIGM settings used are saved in the file `config.txt`.
