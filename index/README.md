Curation of Perturbation-Response Experiments
---------------------------------------------

Our method is dependent on a sufficiently large number publicly available
perturbation experiments that activate or inhibit one of the pathways we were
looking at. The following conditions needed to be met in order for us to
consider an experiment: (1) the compound or factor used for perturbation was
one of our curated list of pathway-perturbing agents (for a list, see table
S0); (2) the perturbation lasted for less than 24 hours to capture
genes that belong to the primary response; (3) there was raw data available for
at least two control arrays and one perturbed array; (4) it was a
single-channel array; (5) we could process the arrays using available
BioConductor packages; (6) the array was not custom-made so we could use
standard annotations.

We curated a list of known pathway activators and inhibitors for 11 pathways,
where the interaction between each compound and pathway is well established
in literature. We then used those as query terms for public perturbation
experiments in the ArrayExpress database 27 and included a total of 223
submissions and 573 experiments in our data set, where each experiment is a
distinct comparison between basal and perturbed arrays. If there were
multiple time points, different cells, different concentrations, or
different perturbing agents within a single database submission, they were
considered as different experiments.
