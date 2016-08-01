Building a Linear Model of Pathway-Response Genes
-------------------------------------------------

For each HGNC symbol, we calculated a model based on mean and standard
deviation of the gene expression level, and computed the z-score as average
number of standard deviations that the expression level in the perturbed array
was shifted from the basal arrays. We then performed LOESS smoothing for all
z-scores in a given experiment using our null model, as described previously
12.

From the z-scores of all experiments and all pathways, we performed a multiple
linear regression with the pathway as input and the z-scores as response
variable for each gene separately:

Z<sub>g</sub> ~ M ... g

Where Zg is the z-score for a given gene g across all input experiments (as a
column vector of experiments). M is a coefficients matrix (rows are
experiments, columns pathways, Fig. 1b) that has the coefficient 1 if the the
experiment had a pathway activated, -1 if inhibited, and 0 otherwise. For
instance, the Hypoxia pathway had experiments with low oxygen conditions set as
1, HIF1A knockdown as -1, and all other experiments as 0. The same is true for
EGFR and EGF treatment vs. EGFR inhibitors respectively, with the additional
coefficients of MAPK and PI3K pathways set to 1 because of known cross-talk
(for a full structure of the cross-talk modeled, see Fig. 1c). As these are
fold changes, we do not allow an intercept.

From the result of the linear model, we selected the top 100 genes per pathway
according to their p-value and took their estimate (the fitted z-scores) as
coefficient. We set all other gene coefficients to 0, so this yielded a matrix
with HGNC symbols in rows and pathways in columns, where each pathway had 100
non-zero gene coefficients (Supplementary Table 9).
