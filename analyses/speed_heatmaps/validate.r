library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')

# load all in ../scores/speed/*.RData
objs = io$load_regex("(.*)\\.RData", "../../scores/speed_nocv")
index = objs$speed_linear$index # index is the same on all objs
# scale the pathway scores for each sample to mean=0, sd=1
sfun = function(o)
    ar$map(o$scores, along=2, function(x) as.numeric(x)) # scale here??
scores = lapply(objs, sfun) %>% ar$intersect_list(along=1)
stopifnot(index$id == rownames(scores[[1]]))

by_method = lapply(scores, function(s) ar$split(s, along=2)) %>%
    b$list$transpose()

# make scores a list of pathways, with matrix cols = methods
pathway2assocs = function(path, sign="activating") {
    path_active = as.numeric(index$pathway == path &
                             index$effect == sign)

    # we use some EGFR exps as positive controls for MAPK/PI3K
    # activity because otherwise we don't have enough activating exps
    map_pi_egfr = c("EGFR.E-GEOD-13168.1",
                    "EGFR.E-GEOD-32975.1",
                    "EGFR.E-GEOD-32975.2")
    if (path %in% c("MAPK","PI3K") && sign=="activating") {
        path_active[index$pathway == "EGFR"] = NA
        path_active[index$id %in% map_pi_egfr] = 1
    }

    score = by_method[[path]]
    st$lm(score ~ path_active) %>%
        arrange(p.value) %>%
        mutate(pathway = path) %>%
        transmute(pathway, sign=sign, method=score,
                  estimate, p.value, size=sum(path_active, na.rm=TRUE))
}
