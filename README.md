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

Recall of perturbations ([figure](figure/2_recall.svg))
-------------------------------------------------------


Functional impact of driver mutations ([figure](figure/3_response_expression.svg))
----------------------------------------------------------------------------------

[analyses/tcga_pathway_per_mutation](analyses/tcga_pathway_per_mutation)

Explaining drug sensitivity ([figure](figure/4_drug_response.svg))
------------------------------------------------------------------

[analyses/drug_assocs](analyses/drug_assocs)

Effect on patient survival ([figure](figure/5_patient_survival.svg))
--------------------------------------------------------------------

[analyses/tcga_survival](analyses/tcga_survival)

![1_overview](https://drive.google.com/uc?id=0B3ZCahnIbtHtWUdFNUQ4bHFPc0E)
![2_recall](https://drive.google.com/uc?id=0B3ZCahnIbtHtTW9XS3IyWjlETTg)
![3_response_expression](https://drive.google.com/uc?id=0B3ZCahnIbtHtd2VQU1k2OUY0Qm8)
![4_drug_response](https://drive.google.com/uc?id=0B3ZCahnIbtHtQ2ttWHh0c2g0cW8)
![5_patient_survival](https://drive.google.com/uc?id=0B3ZCahnIbtHtMHdXUEJWdzZRVGM)
