Microarray Processing
---------------------

Started from the curated list of perturbation-induced gene expression
experiments, we included all single-channel microarrays with at least
duplicates in the basal condition with raw data available that could be
processed by either the limma 28, oligo 29, or affy 30 BioConductor packages
and for which there was a respective annotation package available. Multiple
concentrations or time points in a series of arrays were considered as
individual experiments.

We first calculated a probe-level for 573 full series of arrays, where we
performed quality control of the raw data using RLE and NUSE cutoffs under 0.1
and kept all arrays below that threshold. If after filtering less than two
basal condition arrays remained, the whole experiment was discarded. For the
remaining 568 series we normalized using the RMA algorithm and mapped the probe
identifiers to HGNC symbols.
