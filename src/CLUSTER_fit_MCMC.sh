#!/bin/bash -l

# Standard output and error:
#SBATCH -o results/temp/job-%A.out # Standard output, %A = job ID, %a = job array index
#SBATCH -e results/temp/job-%A.err # Standard error, %A = job ID, %a = job array index

# Job Name:
#SBATCH -J tm

#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=kramer@mpiib-berlin.mpg.de

# --- resource specification (which resources for how long) ---
#SBATCH --partition=general 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000 # memory in MB required by the job
#SBATCH --time=20:00:00 # run time in h:m:s, up to 24h possible

export VIRUS1=$VIRUS1

# --- start from a clean state and load necessary environment modules ---
module purge
module load R/4.0.2

# --- run your executable via srun ---
R --no-save --no-restore <src/fit_MCMC.R >results/Rout/R-tm-$VIRUS1.Rout