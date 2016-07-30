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

The corresponding article for this project is available on bioRxiv:

```
@article {Schubert-PRGs,
	author = {Schubert, Michael and Klinger, Bertram and Kl{\"u}nemann, Martina and Garnett, Mathew J and Bl{\"u}thgen, Nils and Saez-Rodriguez, Julio},
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
[YAML](https://en.wikipedia.org/wiki/YAML) format.

Search terms are supplied in files with name `query.txt`.

Z-scores from gene expression data
----------------------------------

Scripts to download and transform gene expression data, and to generate
z-scores for the perturbation experiments are in the [data](data) directory.

Building the model
------------------

The models we built are available in the [model](model) directory. The one we
used in the publication is called `speed_matrix`.

Computing pathway scores
------------------------

#### On the perturbation experiments

#### On primary tumors of TCGA (The Cancer Gnome Atlas)

TCGA data from [firehose.org](http://firebrowse.org/)

#### On cell lines of the GDSC (Genomics of Drug Sensitivity in Cancer)

GDSC data from [http://www.cancerrxgene.org/gdsc1000](cancerrxgene.org/gdsc1000)

Functional impact of driver mutations
-------------------------------------

[analyses/tcga_pathway_per_mutation](analyses/tcga_pathway_per_mutation)

Explaining drug sensitivity
---------------------------

[analyses/drug_assocs](analyses/drug_assocs)

Effect on patient survival
--------------------------

[analyses/tcga_survival](analyses/tcga_survival)
