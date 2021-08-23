#!/bin/bash

#SBATCH --job-name=random_ranges # the name for the cluster scheduler
#SBATCH --time=2:00:00 # Maximum allowed runtime per iteration
#SBATCH --array=1-1000 # the number of iterations
#SBATCH --mem-per-cpu=8G   # memory per cpu-core
#SBATCH --output=arcus_outputs/rand_ranges%A_%a.out # the name of the output files
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.hemprich-bennett@zoo.ox.ac.uk

module load python/anaconda3/2019.03

source activate $DATA/conda_envs/bat_diet

Rscript array_site_motu95_confidence_intervals.R ${SLURM_ARRAY_TASK_ID}
