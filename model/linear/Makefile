RANDOM_TRIALS = 1000
SETS = pathways \
	   pathways-adrivers \
	   pathways-mut #\
#	   pathways-tdrivers \

FSETS = $(addsuffix .RData, $(SETS))
MOD = $(addprefix mod_, $(FSETS))
BG = $(addprefix mbg_, $(FSETS))
INDEX = $(addprefix index_, $(FSETS))

all: $(MOD) $(BG) $(INDEX)

mod_%.RData: tissue_EN.r
	Rscript $^ --randomize=F --num=1 --sets=$* --outfile=$@

mbg_%.RData: tissue_EN.r
	Rscript $^ --randomize=T --num=$(RANDOM_TRIALS) --sets=$* --outfile=$@

index_%.RData: tissue_index.r mod_%.RData mbg_%.RData
	Rscript $^ $(word 2,$<) $(word 3,$<) $@


