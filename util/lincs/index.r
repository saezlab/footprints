library(dplyr)
b = import('base')
io = import('io')
df = import('data_frame')
lincs = import('data/lincs')

args = commandArgs(TRUE)
INFILES = args[-length(args)] %or% list.files(".", "\\.txt")
OUTFILE = rev(args)[1] %or% "index.RData"

# index should have
#  - control experiments (DMSO treatment)
#  - activity-mod experiments (trt_cp, trt_lig)
#  - rna-mod experiments (trt_oe, trt_sh)
#  - both signs for both
#  - same number of experiments for all (more for control)

control = data.frame(pert_desc="DMSO", pathway="control",
        sign="0", pert_type=c("ctl_vector","ctl_vehicle"))

# read .txt files
pathways = lapply(INFILES, function(fname) {
        message(fname)
        re = io$read_table(fname)
        colnames(re) = c("id", "pert_desc", "pert_type", "sign")
        re$pathway = sub("\\.txt", "", fname)
        re
    }) %>%
    bind_rows() %>%
    bind_rows(control) %>%
    filter(sign %in% c("+","-","0")) %>%
    select(pathway, pert_desc, pert_type, sign)

type_lookup = c(
    ctl_vehicle = "activity",
    ctl_vector = "expression",
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
    filter(!is.na(pert_type), pert_time >= 4, pert_time <= 24)

stopifnot(sum(duplicated(index$distil_id)) == 0)

# save index object
save(index, file="index.RData")
