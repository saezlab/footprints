Pathway-response genes per specific variant
===========================================

#### Rationale

Many variants that occur at a certain frequency in a given cancer type are not
appearing at the same frequency in others. With the different variants occuring
in different cancer types, this can be because they activate differnt pathways
and are selected for or against in either all tissues or a specific one.

#### Goal

Figure out whether common variants have a specific association with a pathway.
This could explain why a variant is selected for in a tissue but not another.

#### Approach

For our set of cancer genes, we select the most common variants across 78 cancer genes (with 20 occcurrences or more) and check if they have a different effect on the pathway-response expression than the rest of the variants.

#### Files

* *gene_variants.r* - assemble the variants of 78 cancer genes using cbioportal
* *variant_impact.r* - use the assembled variants to check which ones are different to the rest
