Pathway-response genes per mutation
===================================

#### Rationale

Many mutations that occur at a certain frequency in a given cancer type are not
appearing at the same frequency in others. 

#### Goal

Here we aim to figure out which mutations give rise to the expression of which
pathway-response genes.

#### Approach

We compare the pathway scores of primary or recurrent tumors that have a
mutation in a given gene with those that have the wild-type version. We only
use information whether a mutation is present or absent, not what it is
exactly. With driver and cna analyses, we exclude all genes that are mutated
but have no predicted functional impact from the analysis.

#### Files

* *snp_drivers.r*, *mutations_annotated{..}.txt* - use SNVs of known cancer
drivers that are linked to at least one pathway
* *snp_all.r* - do not rely on cancer driver list, but use all mutations that
occur above a certain frequency
* *cna_gistic.r*, *cna.txt* - use CNAs instead of SNVs as input; only use
homozygous alterations (this may be subject to change)
