# Use `print-VARNAME` to print arbitrary variables
print-%: ; @echo $* = $($*)

# Always use "all" as default target
.DEFAULT_GOAL := all

# Delete targets on nonzero exit status
.DELETE_ON_ERROR:

# Create rules to help delegate recipes to other Makefiles
#
# There are two parts to this:
# (1) An empty, phony .FORCE target, that will cause
#     external targets to always be built, so that the
#     Makefile there will handle dependencies
# (2) A function that creates explicit rules for each
#     external dependency
#
# Syntax to use this is:
#   $(call ext_dep,DIR_OF_MAKEFILE,TARGETS)
#
# And it will call:
#   make -C DIR_OF_MAKEFILE $(TARGETS)

.PHONY: .FORCE
.FORCE: ;

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
bsub_email = bsub -K -M $(1) -R "rusage[mem=$(1)]"
bsub_log = $(call bsub_email,$(1)) -oo $(addsuffix .log, $(basename $@))
bsub_nolog = $(call bsub_email,$(1)) -oo /dev/null
bsub = bsub_nolog

# Filter strings by a substring they contain
contains = $(foreach v,$(2),$(if $(findstring $(1),$(v)),$(v),))
contains_not = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))
