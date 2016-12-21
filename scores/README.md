Signaling Footprint scores
--------------------------

Each column in the matrix of perturbation-response genes corresponds to a plane
in gene expression space, in which each cell line or tumor sample is located.
If you follow its normal vector from the origin, the distance it spans
corresponds to the pathway score P each sample is assigned (matrix of samples
in rows, pathways in columns). In practice, this is achieved by a simple matrix
multiplication between the gene expression matrix (samples in rows, genes in
columns, values are expression levels) and the model matrix (genes in rows,
pathways in columns, values are our top 100 coefficients):

``` P = E * G ```

We then scaled each pathway or gene set score to have a mean of zero and
standard deviation of one, in order to factor out the difference in strength of
gene expression signatures and thus be able to compare the relative scores
across pathways and samples at the same time.

Pathway and Gene Ontology scores
--------------------------------

We matched our defined set of pathways to the publicly available pathway
databases Reactome 16 and KEGG 31, and Gene Ontology (GO) 17 categories
(Supplementary Tables 1-2), to obtain a uniform set across pathway resources
that makes them comparable. We calculated pathway scores as Gene Set Variation
Analysis (GSVA) scores that are able to assign a score to each individual
sample (unlike GSEA that compares groups).

SPIA scores
-----------

Signaling Pathway Impact Analysis (SPIA) 3 is a method that utilizes the
directionality and signs in a KEGG pathway graph to determine if in a given
pathway structure the available species are more or less available to transduce
a signal. As the species considered for a pathway are usually mRNAs of genes,
this method infers signaling activity by the proxy of gene expression. In order
to do this, SPIA scores require the comparison with a normal condition in order
to compute both their scores and their significance.

We used the SPIA Bioconductor package 3 in our analyses, focussing on a subset
of pathways (Supplementary Table 3). We calculated our scores either for each
cell line compared to the rest of a given tissue where no normals are available
(i.e. for the GDSC and drug response data) or compared to the tissue-matched
normals (for the TCGA data used in driver and survival associations).

Pathifier scores
----------------

As Pathifier 5 requires the comparison with a baseline condition in order to
compute scores, we computed the GDSC/TCGA scores as with SPIA. As gene sets, we
selected Reactome pathways that corresponded to our set of pathways
(Supplementary Table 2), where Pathifier calculated the “signal flow”
from the baseline and compared it to each sample.

PARADIGM scores
---------------

We used the PARADIGM software from the public software repository
(https://github.com/sbenz/Paradigm) and a model of the cell signaling network
32 from the corresponding TCGA publication
(https://tcga-data.nci.nih.gov/docs/publications/coadread_2012/). We normalized
our gene expression data from both GDSC and TCGA using ranks to assign equally
spaced values between 0 and 1 for each sample within a given tissue. We then
ran PARADIGM inference using the same options as in the above publication for
each sample separately. We used nodes in the network representing pathway
activity to our set of pathways (Supplementary Table 4) to obtain pathway
scores that are comparable to the other methods, averaging scores where there
were more than one for a given sample and node.
