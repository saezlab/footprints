# read all .txt files
library(dplyr)
b = import('base')
io = import('io')
df = import('data_frame')
lincs = import('data/lincs')

# read .txt files
pathways = list.files(".", "\\.txt") %>% lapply(function(fname) {
    re = io$read_table(fname)
    colnames(re) = c("id", "pert_desc", "pert_type", "sign")
    re = bind_rows(re, data.frame(id=0, pert_desc="DMSO", pert_type="DMSO", sign="0"))
    re$pathway = sub("\\.txt", "", fname)
    re
}) %>%
    bind_rows() %>%
    filter(sign %in% c("+","0")) %>% #TODO: only pathway activation+controls for now
    select(pathway, pert_desc, sign)

# read CMap inst.info to get experiment ids where this fits
index = lincs$get_index() %>%
    select(distil_id, pert_id, pert_desc, cell_id, pert_time, pert_time_unit, pert_dose, pert_type) %>%
    df$subset(pathways, add_cols=TRUE) %>%
    distinct() %>% #FIXME:
    group_by(pathway, sign) %>%
    do(sample_n(., 100)) %>%
    ungroup() %>%

index = index[!duplicated(index$distil_id),]

# save index object
save(index, file="index.RData")
