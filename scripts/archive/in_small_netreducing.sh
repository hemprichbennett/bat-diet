#!/bin/sh
#$ -cwd     # Use current working directory
#$ -V     # Verbose
#$ -j y     # Maximum output, inc errors
#$ -r y     # Condense error files
#$ -pe smp 1     # Request CPU cores
#$ -l h_rt=150:0:0     # Request runtime (up to 240 hours)
#$ -l h_vmem=10G     # Request RAM per core
#$ -m bea     # Status emails
#$ -t 1-14


module load R



Rscript /data/scratch/btw863/bat-diet/scripts/r/small_netreducing.R ${SGE_TASK_ID}

#Made a smaller version (50 iterations rather than 100) as 100 was timing out for some metrics