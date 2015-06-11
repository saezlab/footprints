# TCGA antibodies with quality scores:
# https://tcga-data.nci.nih.gov/docs/publications/stad_2014/S7.4%20data%20table.pdf
b = import('base')
ar = import('array')

# get icgc RPPA data
icgc = import('data/icgc')
raw = t(icgc$rppa(map.ids = "icgc_specimen_id"))

# discard where many values missing, impute rest
raw = raw[,colSums(is.na(raw)) < 0.8 * nrow(raw)]
raw = raw[rowSums(is.na(raw)) < 0.5 * ncol(raw),]
raw = impute::impute.knn(raw)$data

# map on the pathways, use ratios here
# needed to remove quite a lot, most proteins measurements are sparse
path = list(
    'EGFR' = (raw[,'EGFR_pY1173-R-V'] + raw[,'EGFR_pY992-R-V'] + raw[,'EGFR_pY1068-R-V'])  /
             raw[,'EGFR-R-V'] +
             (raw[,'EGFR_pY1173-R-C'] + raw[,'EGFR_pY1068-R-C']) /
             raw[,'EGFR-R-C'],
    'VEGF' = raw[,'VEGFR2-R-C'] + raw[,'VEGFR2-R-V'],
    'MAPK' = raw[,'MEK1_pS217_S221-R-V'] / raw[,'MEK1-R-V'],
    'JAK-STAT' = raw[,'STAT3_pY705-R-V'],
    'p53' = raw[,'Chk1_pS345-R-C'] / raw[,'Chk1-R-C'],
    'Trail' = raw[,'Bad_pS112-R-V'], # Bcl2, Bcl-xL only prot expr, not comparable
    'NFkB' = raw[,'NF-kB-p65_pS536-R-C'],
    'PI3K' = (raw[,'Akt_pT308-R-V'] + raw[,'Akt_pS473-R-V']) / raw[,'Akt-R-V']
)

# put into matrix, remove outliers
rppa = do.call(cbind, path)
rppa = ar$map(rppa, along=1, function(x) {
    x[x == Inf | x == -Inf] = NA
    x[x > median(x, na.rm=TRUE) + 3*sd(x, na.rm=TRUE) |
      x < median(x, na.rm=TRUE) - 3*sd(x, na.rm=TRUE)] = NA
    scale(x)
})

# save pathway matrix
save(rppa, file="rppa.RData")
