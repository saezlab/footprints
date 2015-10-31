Pathway-response genes per TCGA tissue
======================================

#### Rationale

After we have created the model for pathway-response genes, we calculated the
scores for each pathway and compared those we obtained from tissue-normals with
expression for tumors. This way, we can make statements about increased
expression of pathway-response genes per tumor type (and thus pathway
activation) that does not suffer from generally different expression of those
genes in different tissues.

#### Goal

Using this approach we can see how similar pathway-responsive genes are
activated between different tissues, whether there are sets that are always up-
or down-regulated, and which pathway is important for which cancer type.

#### Approach

We calculated the pathway scores for each cancer type where we have
tissue-matched normal und tumor gene expression available. We calculate scores
for both, and then plot how many standard deviations the tumor score is above
or below the score for the normals.

#### Files

* *tissue_act.r* - calculates pathway scores for each TCGA tissue
* *tact_merge.r* - uses expression data merged between TCGA and GDSC as starting point (to show they merging doesn't change overall scores)
