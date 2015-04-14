MAKEFILEDIR = $(dir $(firstword $(MAKEFILE_LIST)))

PLATFORMS = GPL570 GPL571 GPL96 GPL6244
NORMALISATIONS = rwa none

RWA = $(PLATFORMS:%=frma_rwa_%.ro)
NONE = $(PLATFORMS:%=frma_none_%.ro)

all: $(RWA) $(NONE)

frma_rwa_%.ro: INDEX/index.txt
	bsub -K -R "select[gpfs]" -o /dev/null -M 50000 -R "rusage[mem=50000]" "cd $(MAKEFILEDIR) && \
				 ./normalizeFRMA.r --indexfile=$< \
								   --outfile=$@ \
								   --platform=$* \
								   --normalize=quantile \
								   --summarize=robust_weighted_average"

frma_none_%.ro: INDEX/index.txt
	bsub -K -R "select[gpfs]" -o /dev/null -M 50000 -R "rusage[mem=50000]" "cd $(MAKEFILEDIR) && \
				 ./normalizeFRMA.r --indexfile=$< \
								   --outfile=$@ \
								   --platform=$* \
								   --normalize=none \
								   --summarize=robust_weighted_average"

