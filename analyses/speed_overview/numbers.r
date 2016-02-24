io = import('io')

# speed1 stats: SPEED_db/simple_stats.py in https://github.com/saezlab/speed

records = io$load('../../data/expr.RData')$records
index = io$load('../../data/zscores.RData')$index # the ones we actually use (pass qc, etc)
gatza2014 = io$read_table('../../util/genesets/ng.3073-S2_info.txt', header=TRUE)

speed1 = list(
    pathways = 11,
    experiments = 69,
    contrasts = 215,
    arrays = 572
)

speed2 = list(
    pathways = 19,
    experiments = 202,
    contrasts = 653,
    arrays = 1940
)

speed_matrix = lapply(list(
    pathways = unique(index$pathway),
    experiments = unique(index$accession),
    contrasts = index$id,
    arrays = unique(unlist(lapply(records, function(x) c(x$control, x$perturbed))))
), length)

gatza = list(
    pathways = nrow(gatza2014),
    experiments = length(unique(unlist(stringr::str_match_all(gatza2014$References, "[0-9]+")))),
    contrasts = nrow(gatza2014),
    arrays = NA # uh, manually?
)

df = cbind(speed1, speed2, speed_matrix, gatza)
io$write_table(df, file="numbers.txt", sep="\t")
