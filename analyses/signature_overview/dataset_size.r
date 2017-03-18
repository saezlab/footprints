library(dplyr)
b = import('base')
io = import('io')

# speed1 stats: SPEED_db/simple_stats.py in https://github.com/saezlab/speed

records = io$load('../../data/expr.RData')$records
index = io$load('../../data/zscores.RData')$index # the ones we actually use (pass qc, etc)
gatza2009 = io$read_table('../../util/genesets/ng.3073-S2_info.txt', header=TRUE)
all_records = module_file('../../index') %>%
    list.files("\\.yaml$", recursive=TRUE, full.names=TRUE) %>%
    lapply(function(y) io$read_yaml(y, drop=FALSE)) %>%
    unlist(recursive=FALSE) %>%
    b$omit$null() %>%
    setNames(., lapply(., function(x) x$id))
val_records = module_file('../signature_validation/validation.yaml') %>%
    io$read_yaml()

gsea_speed2016 = list(
    Pathways = 11,
    Datasets = 69,
    Experiments = 215,
    Arrays = 572
)

#speed2 = list(
#    Pathways = 19,
#    Datasets = 202,
#    Experiments = 653,
#    Arrays = 1940
#)

curated = lapply(list(
    Pathways = unique(lapply(all_records, function(x) x$pathway)),
    Datasets = unique(lapply(all_records, function(x) x$accession)),
    Experiments = names(all_records),
    Arrays = all_records %>%
        lapply(function(x) c(x$control, x$perturbed)) %>%
        unlist() %>%
        unique()
), length)

validation = lapply(list(
    Pathways = unique(names(val_records)),
    Datasets = unique(lapply(val_records, function(x) x$expr$accession)),
    Experiments = unique(lapply(val_records, function(x) x$expr$accession)),
    Arrays = val_records %>%
        lapply(function(x) c(x$expr$control, x$expr$perturbed)) %>%
        unlist() %>%
        unique()
), length)

speed_matrix = lapply(list(
    Pathways = unique(index$pathway),
    Datasets = unique(index$accession),
    Experiments = index$id,
    Arrays = records[names(records) %in% index$id] %>%
        lapply(function(x) c(x$control, x$perturbed)) %>%
        unlist() %>%
        unique()
), length)

gsva_gatza = list(
    Pathways = 18, #nrow(gatza2009), # too many duplicates/non-signaling -> use 2009 exps
    Datasets = 1, #length(unique(unlist(stringr::str_match_all(gatza2009$References, "[0-9]+")))),
    Experiments = 18, # nrow(gatza2009),
    Arrays = 287 # by summing up numbers in supp fig. 4
)

df = cbind(measure=names(gsea_speed2016),
           gsea_speed2016,
           gsva_gatza,
           curated,
           speed_matrix,
           validation)

io$write_table(df, file="dataset_size.txt", sep="\t")
