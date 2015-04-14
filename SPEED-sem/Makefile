all: sem_genewise_result.RData

sem_genewise_result.RData: sem_genewise.r
	bsub -K -M 10240 -oo sem_genewise.log -R "rusage[mem=10240]" -R "select[gpfs]" \
		Rscript sem_genewise.r
