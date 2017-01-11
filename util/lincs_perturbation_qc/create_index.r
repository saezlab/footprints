library(dplyr)
b = import('base')
io = import('io')
df = import('data_frame')
lincs = import('data/lincs')

# index should have
#  - control experiments (DMSO treatment)
#  - activity-mod experiments (trt_cp, trt_lig)
#  - rna-mod experiments (trt_oe, trt_sh)
#  - both signs for both
#  - same number of experiments for all (more for control)

# read .txt files
pathways = list.files(".", "\\.txt") %>% lapply(function(fname) {
    re = io$read_table(fname)
    colnames(re) = c("id", "pert_desc", "pert_type", "sign")
    re$pathway = sub("\\.txt", "", fname)
    re
}) %>%
    bind_rows() %>%
    bind_rows(data.frame(pert_desc="DMSO", pathway="control", sign="0")) %>%
    filter(sign %in% c("+","-","0")) %>%
    select(pathway, pert_desc, sign)

type_lookup = c(
    ctl_vehicle = "control",
    trt_cp = "activity",
    trt_lig = "activity",
    trt_oe = "expression",
    trt_sh = "expression"
)

# read CMap inst.info to get experiment ids where this fits
index = lincs$get_index() %>%
    select(distil_id, pert_id, pert_desc, cell_id, pert_time, pert_time_unit, pert_dose, pert_type) %>%
    filter(pert_time_unit == "h") %>%
    df$subset(pathways, add_cols=TRUE) %>%
    mutate(pert_type = type_lookup[pert_type]) %>%
    filter(!is.na(pert_type)) %>%
#    group_by(pathway, sign, pert_type) %>%
#    do(sample_n(., 100, replace=TRUE)) %>% #FIXME: no replace, but let s<n be?
#    ungroup()

# save index object
save(index, file="index.RData")
