# get icgc RPPA data
# save in same format as other scores
icgc = import('data/icgc')
rppa = t(icgc$rppa(map.ids = "icgc_specimen_id"))
save(rppa, file="rppa.RData")

#TODO: mapping on the pathways here
