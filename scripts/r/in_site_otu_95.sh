#!/bin/sh
#$ -cwd     # Use current working directory
#$ -V     # Verbose
#$ -j y     # Maximum output, inc errors
#$ -r y     # Condense error files
#$ -pe smp 1     # Request CPU cores
#$ -l h_rt=100:0:0     # Request runtime (up to 240 hours)
#$ -l h_vmem=5G     # Request RAM per core
#$ -m bea     # Status emails
#$ -t 1-11


module load R



Rscript /data/scratch/btw863/bat-diet/scripts/r/array_site_motu95_confidence_intervals.R ${SGE_TASK_ID}
