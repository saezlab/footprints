Signaling footprints via pathway-responsive genes (PRGs)
========================================================

Numerous pathway methods have been developed to quantify the signaling state of
a cell from gene expression data, usually from the abundance of transcripts of
pathway members, and are hence unable to take into account post-translational
control of signal transduction. Gene expression signatures of pathway
perturbations can capture this, but they are closely tied to the experimental
conditions that they were derived from. We overcome both limitations by
leveraging a large compendium of publicly available perturbation experiments to
define consensus signatures for pathway activity. We find that although
individual expression signatures are heterogeneous, there is a common core of
responsive genes that describe pathway activation in a wide range of
conditions. These signaling footprints better recover pathway activity than
existing methods and provide more meaningful associations with (i) known driver
mutations in primary tumors, (ii) drug response in cell lines, and (iii)
survival in cancer patients, making them more suitable to assess the activity
status of signaling pathways.

The corresponding article for this project is available on
[bioRxiv](http://biorxiv.org/content/early/2016/07/25/065672)
([pdf](http://biorxiv.org/content/early/2016/07/25/065672.full.pdf)).

```
@article {Schubert-PRGs,
	author = {Schubert, Michael and Klinger, Bertram and Kl{\"u}nemann, Martina and 
              Garnett, Mathew J and Bl{\"u}thgen, Nils and Saez-Rodriguez, Julio},
	title = {Perturbation-response genes reveal signaling footprints in cancer gene expression},
	year = {2016},
	doi = {10.1101/065672},
	publisher = {Cold Spring Harbor Labs Journals},
	URL = {http://biorxiv.org/content/early/2016/07/25/065672},
	eprint = {http://biorxiv.org/content/early/2016/07/25/065672.full.pdf},
	journal = {bioRxiv}
}
```

Perturbation experiments
------------------------

The 11 pathways comprised of 208 submissions to
[ArrayExpress](https://www.ebi.ac.uk/arrayexpress/) with a total of 580
experiments are available in the [index](index) directory in
[YAML](https://en.wikipedia.org/wiki/YAML) format. We consider the following
pathways:

 * EGFR
 * Hypoxia
 * JAK-STAT
 * MAPK
 * NFkB
 * p53/DDR (DNA damage response)
 * PI3K
 * TGFb
 * Trail
 * VEGF (and PDGF)

Search terms are supplied in files with name `query.txt`. Files and experiments
we excluded (because of QC failure or we were not sure if the perturbation
corresponds to the phenotype we want to observe) are indicated using the
`.excluded` suffix or a commented entry in the file. Reason for exclusion is
mentioned in the files.

Z-scores from gene expression data
----------------------------------

Scripts to download and transform gene expression data, and to generate
z-scores for the perturbation experiments are in the [data](data) directory.

We normalized and QC'd each series as whole (`normalize_data.r`), then
assembled the relevant expression data (`expr.r`) and computed z-scores for
each perturbation experiment (`zscores.r`).

Building the model
------------------

The models we built are available in the [model](model) directory. The one we
used in the publication is called `speed_matrix`.

For this, we fit a linear model on the z-scores with a binary matrix indicating
pathway perturbations (incl. perturbations of multiple pathways) as the
independent variable. We select the 100 most significant genes and use their
z-scores as coefficients in the model.

Computing pathway scores
------------------------

#### Different pathway methods considered

We computed [pathway scores](scores) for the PRGs, and the [corresponding
pathways](config/pathway_mapping.yaml) for Gene Ontology and Reactome genesets
(using GSVA), Signaling Pathway Impact Analysis (SPIA;
[article](http://bioinformatics.oxfordjournals.org/content/25/1/75.short)
and [R package](http://bioconductor.org/packages/release/bioc/html/SPIA.html)),
Pathifier ([article](http://www.pnas.org/content/110/16/6388.short),
[R package](http://bioconductor.org/packages/release/bioc/html/pathifier.html)),
and PARADIGM
([article](http://bioinformatics.oxfordjournals.org/content/26/12/i237.short),
[tool](https://github.com/sbenz/Paradigm) using the [TCGA signaling
network](https://tcga-data.nci.nih.gov/docs/publications/coadread_2012/)).

#### On the perturbation experiments

For the [pathway scores derived from perturbations](scores/speed), we used the
fold changes (PRGs, Gene Ontology, Reactome) or the basal samples and perturbed
samples as control and perturbed, respectively (SPIA, Pathifier).

#### On primary tumors of TCGA (The Cancer Genome Atlas)

We took TCGA data from [firehose.org](http://firebrowse.org/) and computed all
pathway scores for primary tumors (`01A` in barcode) where a tissue-matched
normal (`11A`) was available (required for SPIA and Pathifier).

#### On cell lines of the GDSC (Genomics of Drug Sensitivity in Cancer)

We took the GDSC data from
[cancerrxgene.org/gdsc1000](http://www.cancerrxgene.org/gdsc1000)
([article](http://www.sciencedirect.com/science/article/pii/S0092867416307462)),
computing pathway scores for each cell line (for SPIA and Pathifier compared to
all other cell lines of the same TCGA label).

Overview of the approach
------------------------

<img src=https://drive.google.com/uc?id=0B3ZCahnIbtHtWUdFNUQ4bHFPc0E width=50%>

A. Reasoning about pathway activation. Most pathway approaches make use of
either the set (top panel) or network (middle panel) of signaling molecules to
make statements about a possible activation, while our approach considered the
genes affected by perturbing them.
B. Workflow of data curation and model building. (1) Finding and curation of
208 publicly available experiment series in the ArrayExpress database, (2)
    Extracting 556 perturbation experiments from series' raw data, (3)
    Performing QC metrics and discarding failures, (4) Computing z-scores per
    experiment, (5) Using a multiple linear regression to fit genes responsive
    to all pathways simultaneously obtaining the z-coefficients matrix, (6)
    Assigning pathway scores using the coefficients matrix and basal expression
    data. See methods section for details. Image credit Supplementary Note 1.
    C. Structure of the perturbation-response model. For the multiple linear
    regression, we set the coefficients or perturbed pathways to 1 if a pathway
    was perturbed, 0 otherwise. In addition, EGFR perturbation also had MAPK
    and PI3K coefficients set, and TNFa had NFkB set.

Recall of perturbations ([analyses/speed_raw](analyses/speed_raw))
-------------------------------------------------------

<img src=https://drive.google.com/uc?id=0B3ZCahnIbtHtTW9XS3IyWjlETTg width=80%>

**A.** T-SNE plots for separation of perturbation experiments with different
pathway perturbations in different colors. Fold changes of genes in individual
perturbation experiments (10% FDR) do not cluster by pathway (left). Using a
consensus signature of genes whose z-score is most consistently deregulated for
each pathway instead, we can observe distinct clusters of perturbed pathways
(right). Details Supplementary Note 2.
**B.** Associations between perturbed pathways and the scores obtained by the model
of pathway-responsive genes (PRGs). Along the diagonal each pathway is strongly
(p&lt;10<sup>-10</sup>) associated with its own perturbation. Significant off-diagonal
elements are sparse and only occur (p&lt;10<sup>-5</sup>) where there is biologically known
cross-activation.
**C.** Heatmap of relative pathway scores in each perturbation experiment. 523
experiments in columns, annotated with the perturbation effect (green for
activation, orange for inhibition) and pathway perturbed (same order as
b). Pathway scores in rows cluster between EGFR/MAPK and to a
lesser extent PI3K, and TNFa/NFkB. Color indicates activation or
inhibition strength.
**D.** ROC curves for different methods ranking perturbation experiments by
their pathway score. PRGs show better performance for all pathways
except JAK-STAT and NFkB, where other methods are equal. Gene Ontology
and Reactome scores obtained by Gene Set Variation Analysis (GSVA).
Pathifier using Reactome gene sets.
**E.** Correlation of pathway scores in basal gene expression of cell lines
in the GDSC panel. Positive correlation in blue, negative in red.
Circle size and shade correspond to correlation strength. Pathways that
showed cross-activation in point b are more highly correlated in basal
expression as well. 
**F.** Stability of basal pathway scores when bootstrapping input
experiments. The variance of pathway scores in cell lines given
bootstraps more than five times as high compared to the variance of
bootstraps given cell lines for all pathways except two (Trail and
VEGF), where it is roughly twice as high.

Functional impact of driver mutations ([analyses/tcga_pathway_per_mutation](analyses/tcga_pathway_per_mutation))
----------------------------------------------------------------------------------

<img src=https://drive.google.com/uc?id=0B3ZCahnIbtHtd2VQU1k2OUY0Qm8 width=80%>

**A.** Volcano plot of pan-cancer associations between driver mutations and copy
number aberrations with differences in pathway score. Pathway scores calculated
from basal gene expression in the TCGA for primary tumors. Size of points
corresponds to occurrence of aberration. Type of aberration is indicated by
superscript “mut” if mutated and “amp”/”del” if
amplified or deleted, with colors as indicated. Effect sizes on the horizontal
axis larger than zero indicate pathway activation, smaller than zero inferred
inhibition. P-values on the vertical axis FDR-adjusted with a significance
threshold of 5%. Associations shown without correcting for different cancer
types. Associations with a black outer ring are also significant if corrected.
**B.** Comparison of pathway scores (vertical axes) across different methods
(horizontal axes) for TP53 and KRAS mutations, EGFR amplifications and VHL
mutations. Wald statistic shown as shades of green for downregulated and red
for upregulated pathways. P-value labels shown as indicated. White squares
where a pathway was not available for a method.

Explaining drug sensitivity ([analyses/drug_assocs](analyses/drug_assocs))
------------------------------------------------------------------

<img src=https://drive.google.com/uc?id=0B3ZCahnIbtHtQ2ttWHh0c2g0cW8 width=80%>

**A.** Volcano plot of pan-cancer associations between PRG pathway scores and drug
response (log10 IC50). Pathway scores computed using basal gene expression in
GDSC cell lines. Associations corrected for cancer type. Size of points
corresponds to number of cell lines screened with a particular drug. Effect
size corresponds to 10-fold change in IC50 per standard deviation of the
pathway score. Values smaller than zero indicate sensitivity markers (green)
and greater than zero resistance markers (red). P-values FDR-corrected.
**B.** Pathway context of the strongest associations between EGFR/MAPK pathways
and their inhibitors.
**C.** Comparison of the associations obtained by different pathway methods.
Number of associations on the vertical, FDR on the horizontal axis. PRGs
yield more and stronger associations than all other pathway methods.
Mutation associations are only stronger for TP53/Nutlin-3a and drugs that
were specifically designed to bind to a mutated protein. PARADIGM not shown
because no associations < 10% FDR.
**D.** Comparison of stratification by mutations and pathway scores. MAPK
pathway (BRAF, NRAS, or KRAS) mutations and Trametinib on top left panel,
AZ628 bottom left, BRAF mutations and Dabrafenib top right, and p53
pathway/TP53 mutations/Nutlin-3a bottom right. For each of the four cases,
the leftmost violin plot shows the distribution of IC50s across all cell
lines, followed by a stratification in wild-type (green) and mutant cell
lines (blue box). The three rightmost violin plots show stratification of
all the cell lines by the top, the two middle, and the bottom quartile of
inferred pathway score (indicated by shade of color). The two remaining
violin plots in the middle show mutated (BRAF, KRAS, or NRAS; blue color)
or wild-type (TP53; green color) cell lines stratified by the top- and
bottom quartiles of MAPK or p53 pathways scores (Mann-Whitney U test
statistics as indicated).

Effect on patient survival ([analyses/tcga_survival](analyses/tcga_survival))
--------------------------------------------------------------------

<img src=https://drive.google.com/uc?id=0B3ZCahnIbtHtMHdXUEJWdzZRVGM width=70%>

**A.** Pan-cancer associations between pathway scores and patient survival.
Pathways on the horizontal, different methods on the vertical axis.
Associations of survival increase (green) and decrease . Significance labels as
indicated. Shades correspond to effect size, p-values as indicated.
**B.** Volcano plot of cancers-specific associations between patient survival and
inferred pathway score using PRGs. Effect size on the horizontal axis. Below
zero indicates increased survival (green), above decreased survival (red).
FDR-adjusted p-values on the vertical axis. Size of the dots corresponds to
number of patients in each cohort.
**D.** Kaplan-Meier curves of individual associations for kidney (KIRC), low-grade
glioma (LGG) and adrenocortical carcinoma (ACC). Pathway scores are split in
top- and bottom quartiles and center half. Lines show the fraction of patients
(vertical axis) that are alive at a given time (horizontal axis) within one
year. P-values for discretized scores.
