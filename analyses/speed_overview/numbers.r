io = import('io')

# speed1 stats: SPEED_db/simple_stats.py in https://github.com/saezlab/speed

records = io$load('../../data/expr.RData')$records
index = io$load('../../data/zscores.RData')$index # the ones we actually use (pass qc, etc)
gatza2014 = io$read_table('../../util/genesets/ng.3073-S2_info.txt', header=TRUE)

speed1 = list(
    Pathways = 11,
    Datasets = 69,
    Experiments = 215,
    Arrays = 572
)

speed2 = list(
    Pathways = 19,
    Datasets = 202,
    Experiments = 653,
    Arrays = 1940
)

speed_matrix = lapply(list(
    Pathways = unique(index$pathway),
    Datasets = unique(index$accession),
    Experiments = index$id,
    Arrays = unique(unlist(lapply(records, function(x) c(x$control, x$perturbed))))
), length)

gatza = list(
    Pathways = nrow(gatza2014),
    Datasets = length(unique(unlist(stringr::str_match_all(gatza2014$References, "[0-9]+")))),
    Experiments = nrow(gatza2014),
    Arrays = NA # uh, manually?
)

df = cbind(measure=names(speed1), speed1, speed2, speed_matrix, gatza)
io$write_table(df, file="numbers.txt", sep="\t")
