Associations with known driver mutations and CNAs
-------------------------------------------------

For comparing the impact of mutations across different pathway methods, we used
TCGA cohorts where tissue-matched controls were available, leaving 6549 samples
across 13 cancer types. For mutated genes, we considered all genes that had a
change of coding sequence (SNP, small indels in MAF files) as mutated and all
others as not mutated. For copy number alterations (CNAs), we used the
thresholded GISTIC 33 scores, where we considered homozygous deletions (-2) and
strong amplifications (2) as altered, no change (0) as basal and discarded
intermediate values (-1, 1) from our analysis. We focussed our analysis of the
mutations and copy number alterations on the subset of 464 driver genes that
were also used in the GDSC. We used the sets of mutations and CNAs to compute
the linear associations between samples for all different methods we looked at.
We did not regress out the cancer type in order to keep associations where
mutations/CNAs are highly correlated with it, but highlighted all associations
that passed the significance threshold of FDR<5% (for each pathway method
individually) after such a correction.

#### Files

* *snp_drivers.r*, *mutations_annotated{..}.txt* - use SNVs of known cancer
drivers that are linked to at least one pathway
* *snp_all.r* - do not rely on cancer driver list, but use all mutations that
occur above a certain frequency
* *cna_gistic.r*, *cna.txt* - use CNAs instead of SNVs as input; only use
homozygous alterations (this may be subject to change)
