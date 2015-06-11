# TCGA antibodies with quality scores:
# https://tcga-data.nci.nih.gov/docs/publications/stad_2014/S7.4%20data%20table.pdf

# get icgc RPPA data
icgc = import('data/icgc')
raw = t(icgc$rppa(map.ids = "icgc_specimen_id"))

# map on the pathways, use ratios here
path = list(
    'EGFR' = (raw[,'EGFR_pY1173-R-V'] + raw[,'EGFR_pY992-R-V'] + raw[,'EGFR_pY1068-R-V'])  /
             raw[,'EGFR-R-V'] +
             (raw[,'EGFR_pY1173-R-C'] + raw[,'EGFR_pY1068-R-C']) /
             raw[,'EGFR-R-C'] +
             (raw[,'EGFR_pY1173'] + raw[,'EGFR_pY1068']) /
             raw[,'EGFR'],
    'VEGF' = raw[,'VEGFR2'] + raw[,'VEGFR2-R-C'] + raw[,'VEGFR2-R-V'],
    'MAPK' = (raw[,'MEK1_pS217_S221'] + raw[,'MEK1_pS217_S221-R-V']) /
             (raw[,'MEK1'] + raw[,'MEK1-R-V']),
    'JAK-STAT' = raw[,'STAT3_pY705'] + raw[,'STAT3_pY705-R-V'],
    'p53' = (raw[,'Chk1_pS345'] / raw[,'Chk1']) +
            (raw[,'Chk1_pS345-R-C'] + raw[,'Chk1-R-C']),
    'Trail' = raw[,'Bad_pS112'] + raw[,'Bad_pS112-R-V'], # Bcl2, Bcl-xL only prot expr, not comparable
    'NFkB' = raw[,'NF-kB-p65_pS536'] + raw[,'NF-kB-p65_pS536-R-C'],
    'PI3K' = (raw[,'Akt_pT308'] + raw[,'Akt_pS473']) / raw[,'Akt'] +
             (raw[,'Akt_pT308-R-V'] + raw[,'Akt_pS473-R-V']) / raw[,'Akt-R-V']
)
rppa = do.call(cbind, path)

# save pathway matrix
save(rppa, file="rppa.RData")
