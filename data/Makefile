METHODS = none combatgpl combatgse dwdgpl dwdgse
RMA = $(METHODS:%=SPEED2mat_rma_%.RData)
FRMA = $(subst rma,frma,$(RMA))
SOFT = $(subst rma,soft,$(RMA))

all: rma frma

rma: $(RMA)

frma: $(FRMA)

soft: $(SOFT)


SPEED2mat_%_dwdgpl.RData: SPEED2mat.r
	Rscript SPEED2mat.r --norm=$* --method=dwd --batch=gpl

SPEED2mat_%_dwdgse.RData: SPEED2mat.r
	bsub -K -M 20480 -R "select[gpfs]" -R "rusage[mem=20480]" \
		Rscript SPEED2mat.r --norm=$* --method=dwd --batch=gse


SPEED2mat_%_combatgpl.RData: SPEED2mat.r
	Rscript SPEED2mat.r --norm=$* --method=combat --batch=gpl

SPEED2mat_%_combatgse.RData: SPEED2mat.r
	Rscript SPEED2mat.r --norm=$* --method=combat --batch=gse

SPEED2mat_%_none.RData: SPEED2mat.r
	Rscript SPEED2mat.r --norm=$* --method=none


#%.RData: elasticNet.r lists
SPEED2vecs.RData: $(EXPR)
	bsub -K -R "select[gpfs]" -oo make.log.vecs Rscript SPEED2vecs.r $@




soft: SPEED1_soft.RData SPEED2_soft.RData

SPEED%_soft.RData: ../SPEED-API/other_scripts/SPEED%.db
	Rscript processSoft2RData.r --db $^ --outfile $@


