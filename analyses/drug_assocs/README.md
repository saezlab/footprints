Drug associations using GDSC cell lines
---------------------------------------

We performed drug association using an ANOVA between 265 drug IC50s and 11
inferred pathway scores conditioned on MSI status, doing a total of 2915
comparisons for which we correct the p-values using the false discovery rate.
For pan-cancer associations, we used the cancer type as a covariate in order to
discard the effect that different tissues have on the observed drug response.
While this will also remove genuine differences in pathway activation between
different cancer types, we would not be able to distinguish between those and
other confounders that impact the sensitivity of a certain cell line from a
given tissue to a drug. Our pan-cancer association are thus the same of
intra-tissue differences in drug response explained by inferred (our method,
GO, or Reactome) pathway scores.

We also selected two of our strongest associations to investigate whether they
provide additional information of what is known by mutation data. For two MEK
inhibitors, we show the difference between wild-type and mutant MAPK pathway
(defined as a mutation in either MAP2K1, NRAS, KRAS, or BRAF) with a
discretized pathway score (upper and lower quartile vs. the rest), as well as
the combination between the upper quartile of tissue-specific pathway scores
and presence of a MAPK mutation.
