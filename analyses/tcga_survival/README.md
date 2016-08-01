Survival associations using TCGA data
-------------------------------------

Starting from the pathway scores derived with GO/Reactome GSEA, SPIA,
Pathifier, PARADIGM, and our method on the TCGA data as described
above, we used Cox Proportional Hazard model (R package survival) to
calculate survival associations for pan-cancer and each
tissue-specific cohort. For the pan-cancer cohort, we regressed out
the effect of the study and age of the patient, and fitted the more
for each pathway and method used. For the tissue-specific cohorts, we
regressed out the age of the patients. We adjusted the p-values using
the FDR method for each method and for each method and study
separately. We selected a significance threshold of 5 and 10% for the
pan-cancer and cancer-specific associations for which we show a matrix
plot and a volcano plot of associations, respectively.

In order to get distinct classes needed for interpretable Kaplan-Meier survival
curves (Fig. 4c), we split all obtained pathway scores in upper, the two
middle, and lower quartile and respectively to show for the three examples of
associations found.
