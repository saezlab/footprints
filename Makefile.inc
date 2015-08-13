# target word, split by '.', '/' and '_'
tw = $(word $(1), $(subst /, ,$(subst ., ,$(subst _, ,$@))))

# will return every combination of 2 variables separated by '_'
grid = $(foreach x,$(1),$(foreach y,$(2),$(x)_$(y)))

# job submission with resources
bsub = bsub -K -M $(1) -R "rusage[mem=$(1)]" -R "select[panfs_nobackup_research]"
bsub_log = $(call bsub,$(1)) -oo $(@:%.RData=%.log)
