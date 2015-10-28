# Use `print-VARNAME` to print arbitrary variables
print-%: ; @echo $* = $($*)

# Delete targets on nonzero exit status
.DELETE_ON_ERROR:

# Create rules to help delegate recipes to other Makefiles
#
# There are two parts to this:
# (1) An empty, phony .FORCE target, that will cause
#     external targets to always be built, so that the
#     Makefile there will handle dependencies
# (2) A default rule that uses the Makefile in the target
#     directory to build a file for all files where no
#     target is specified.
.PHONY: .FORCE
.FORCE: ;

#.DEFAULT: .FORCE
#	@echo running: $@
#	make -C $(dir $@) $(notdir $@)


# Alternative approach to above
# this will only work if file to be made is same dir as Makefile
define EXT_DEP

$1/$2: .FORCE
	make -C $(1) $(2)

endef

ext_dep = $(foreach tt,$(2),$(eval $(call EXT_DEP,$(1),$(tt))))


# Tokenizer for the target variable: split by '.', '/' and '_'
tw = $(word $(1), $(subst /, ,$(subst ., ,$(subst _, ,$@))))

# Will return every combination of 2 variables separated by '_'
grid = $(foreach x,$(1),$(foreach y,$(2),$(x)_$(y)))

# Job submission with resources
bsub = bsub -K -M $(1) -R "rusage[mem=$(1)]" -R "select[panfs_nobackup_research]"
bsub_log = $(call bsub,$(1)) -oo $(@:%.RData=%.log)
