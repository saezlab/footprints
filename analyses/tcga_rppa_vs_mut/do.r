io = import('io')
ar = import('array')
st = import('stats')
tcga = import('data/tcga')
mut = tcga$mutations()

# summarise the "mapk-activating" mutations
has_mapk = mut %>%
    filter(Hugo_Symbol %in% c("HRAS","KRAS","NRAS","BRAF","CRAF","MAPK2K1")) %>%
    select(Tumor_Sample_Barcode, study) %>%
    distinct() %>%
    mutate(mut = TRUE)
no_mapk = mut %>%
    filter(! Tumor_Sample_Barcode %in% has_mapk$Tumor_Sample_Barcode) %>%
    select(Tumor_Sample_Barcode, study) %>%
    distinct() %>%
    mutate(mut = FALSE)

# total number of mapk mutated/not mutated samples available, disregarding speed/rppa
merge(has_mapk %>% group_by(study) %>% summarise(mut = n()),
      no_mapk %>% group_by(study) %>% summarise(no_mapk = n()),
      by = "study")

# prepare fur subsetting
mapk = rbind(has_mapk, no_mapk)
mapk = mapk[!duplicated(substr(mapk$Tumor_Sample_Barcode, 1, 16)),]
rownames(mapk) = substr(mapk$Tumor_Sample_Barcode, 1, 16)

# get data where we have mutational status, speed scores, and rppa data
speed = io$load("../../scores/tcga/speed_linear.RData")[,"MAPK"]
names(speed) = substr(names(speed), 1, 16)
rppa = t(tcga$rppa())[,"MAP2K1|MEK1_pS217_S221"]
rppa = rppa[!is.na(rppa)]
names(rppa) = substr(names(rppa), 1, 16)
ar$intersect(speed, rppa, mapk)

# how many mutated samples we have for each tissue
merge(mapk %>% group_by(study) %>% summarise(sum(mut)),
      mapk %>% group_by(study) %>% summarise(sum(!mut)),
      by="study")

mapk$speed = speed
mapk$rppa = rppa

mat = data.matrix(mapk[-c(1:3)])
st$lm(mat ~ mapk$mut) %>%
    filter(term != "(Intercept)")

#TODO: other way around, which mutations activate mapk in breast cancer?
#  which mutations are used/planned to be used to stratify patients for MEKi treatment?
